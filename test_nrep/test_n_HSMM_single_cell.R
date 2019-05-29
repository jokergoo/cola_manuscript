library(GetoptLong)

GetoptLong(
	"run=i", "run",
	"rep=i", "rep"
)

options(showWarnCalls = TRUE, showErrorCalls = TRUE)

# root = "/home/guz"
BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"
root = "/desktop-home/guz"


library(cola)
library(GetoptLong)
library(RColorBrewer)

library(HSMMSingleCell)

data(HSMM_expr_matrix)
data(HSMM_sample_sheet)

# `HSMM_expr_matrix` is a FPKM matrix
m = adjust_matrix(log10(HSMM_expr_matrix + 1))
anno = HSMM_sample_sheet[, c("Hours", "Media", "State")]
anno_col = list(
	Hours = structure(brewer.pal(9, "Blues")[c(2, 4, 6, 8)], names = c("0", "24", "48", "72")),
	Media = c("GM" = "orange", "DM" = "purple"),
	State = c("1" = "red", "2" = "blue", "3" = "green"))

gt = readRDS(qq("@{root}/project/development/cola_examples/HSMM_single_cell/gene_type_gencode_v17.rds"))
m = m[gt[rownames(m)] == "protein_coding", , drop = FALSE]

register_NMF()

rl = run_all_consensus_partition_methods(
	m,
	top_n = c(1000, 2000, 3000, 4000), 
	partition_repeat = rep,
	mc.cores = 4, 
	anno = anno,
	anno_col = anno_col
)
saveRDS(rl, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/HSMM_single_cell_subgroup_repeat_@{rep}_run_@{run}.rds"))

q(save = "no")

# for(rep in c(25, 50, 100, 200)) {
# 	for(run in 1:100) {
# 		cmd = qq("module load R/3.3.1; Rscript /icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_n_HSMM_single_cell.R --rep @{rep} --run @{run}")
# 		dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/")
# 		cmd = qq("perl /desktop-home/guz/project/development/ngspipeline2/bsub_single_line.pl --hour 50 --memory 20 --core 4 --name HSMM_single_cell_subgroup_repeat_@{rep}_run_@{run} --dir @{dir} --command '@{cmd}'")
# 		system(cmd)
# 	}
# }


library(GetoptLong)
library(cola)
library(circlize)

for(nrep in c(25, 50, 100, 200)) {

df_all = NULL
for(i in 1:100) {
	qqcat("loading @{i}/100...\n")
	try({
	res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/HSMM_single_cell_subgroup_repeat_@{nrep}_run_@{i}.rds"))

	df_all = rbind(df_all, do.call(rbind, lapply(2:6, function(k) {
		df = get_stats(res, k = k)
		df = as.data.frame(df)
		df$method = rownames(df)
		df$i_simulate = i
		best_k = guess_best_k(res)[rownames(df), ]
		df$best_k = best_k$best_k
		df$best_k_sig = best_k[, ncol(best_k)]
		rownames(df) = NULL
		df$PAC = 1 - df$PAC
		names(df)[names(df) == "PAC"] = "1-PAC"
		df
	})))

	})
}

pdf(qq("@{BASE_DIR}/image/HSMM_single_cell_subgroup_repeat_@{nrep}_mean_vs_vc.pdf"))
par(mfrow = c(2, 2))
for(stat in c("cophcor", "1-PAC", "mean_silhouette", "concordance")) {
	x = tapply(df_all[[stat]], paste0(df_all$method, df_all$k), mean)
	y = tapply(df_all[[stat]], paste0(df_all$method, df_all$k), sd)
	col = tapply(df_all$k, paste0(df_all$method, df_all$k), unique)
	is_sig = tapply(seq_len(nrow(df_all)), paste0(df_all$method, df_all$k), 
		function(i) sum(df_all[i, "1-PAC"] >=0.9)/length(i) > 0.5)

	y = y/x
	x = x^4
	plot(x, y, col = add_transparency(col, ifelse(is_sig, 0, 0.75)), 
		pch = ifelse(is_sig, 1, 4),
		cex = ifelse(is_sig, 1, 0.6),
		main = stat,
		xlim = c(0, 1), axes = FALSE, xlab = "mean", ylab = "variance coefficient")
	axis(side = 1, at = c(0, 0.4, 0.6, 0.8, 0.9, 0.95, 1)^4, labels = c(0, 0.4, 0.6, 0.8, 0.9, 0.95, 1))
	axis(side = 2)
	# points(x["ATC:mclust2"], y["ATC:mclust2"], pch = 16)
	od = order(x)
	fit = loess(y[od] ~ x[od])
	lines(x[od], predict(fit, x[od]), col = "grey")
	box()
}
dev.off()

}

get_vaules = function(stat) {

	m = NULL
	for(nrep in c(25, 50, 100, 200)) {
		df_all = NULL
		for(i in 1:100) {
			qqcat("loading HSMM_single_cell_subgroup_repeat_@{nrep}_run_@{i}.rds...\n")
			res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/HSMM_single_cell_subgroup_repeat_@{nrep}_run_@{i}.rds"))

			df_all = rbind(df_all, do.call(rbind, lapply(2:6, function(k) {
				df = get_stats(res, k = k)
				df = as.data.frame(df)
				df$method = rownames(df)
				df$i_simulate = i
				best_k = guess_best_k(res)[rownames(df), ]
				df$best_k = best_k$best_k
				df$best_k_sig = best_k[, ncol(best_k)]
				rownames(df) = NULL
				df$PAC = 1 - df$PAC
				names(df)[names(df) == "PAC"] = "1-PAC"
				df
			})))
		}
		x = tapply(df_all[[stat]], paste0(df_all$method, df_all$k), mean)
		y = tapply(df_all[[stat]], paste0(df_all$method, df_all$k), sd)
		k = tapply(df_all$k, paste0(df_all$method, df_all$k), unique)
		
		m = rbind(m, cbind(x = x, y = y, k = k, nrep = nrep))
	}
	as.data.frame(m)
}

for(stat in c("1-PAC", "mean_silhouette", "cophcor", "concordance")) {
	m = get_vaules(stat)

	library(ggplot2)
	library(scales)
	pdf(qq("@{BASE_DIR}/image/HSMM_single_cell_subgroup_repeat_@{stat}_mean_vs_vc.pdf"), width = 15, height = 5)
	p = ggplot(m, aes(x = x, y = y/x)) +
		geom_point(aes(col = factor(k))) + geom_smooth(se = FALSE, method = "loess", col = "grey", lwd = 1) +
		facet_wrap(~ factor(nrep, levels = c(25, 50, 100, 200)), nrow = 1) +
		xlab(qq("mean @{stat}")) + ylab("variance coefficient")
	print(p)
	par(mfrow = c(1, 4))
	plot(m[m$nrep == 25, 1], m[m$nrep == 50, 1], col = (m[m$nrep == 25, "k"]), xlab = "nrep = 25", ylab = "nrep = 50", main = qq("mean @{stat}"))
	plot(m[m$nrep == 50, 1], m[m$nrep == 100, 1], col = (m[m$nrep == 50, "k"]), xlab = "nrep = 50", ylab = "nrep = 100", main = qq("mean @{stat}"))
	plot(m[m$nrep == 100, 1], m[m$nrep == 200, 1], col = (m[m$nrep == 100, "k"]), xlab = "nrep = 100", ylab = "nrep = 200", main = qq("mean @{stat}"))
	plot(m[m$nrep == 25, 1], m[m$nrep == 200, 1], col = (m[m$nrep == 25, "k"]), xlab = "nrep = 25", ylab = "nrep = 200", main = qq("mean @{stat}"))
	dev.off()
}

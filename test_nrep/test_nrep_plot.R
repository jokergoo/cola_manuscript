

library(GetoptLong)
library(cola)
library(circlize)

GetoptLong(
	"dataset=s", "dataset"
)

# root = "/home/guz"
BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"

# get the mean/sd of stat in each combination of method-k	
get_values = function(stat) {
	m = NULL
	df_return = NULL
	for(nrep in c(25, 50, 100, 200)) {
		df_all = NULL
		for(i in 1:100) {
			qqcat("loading @{dataset}_subgroup_repeat_@{nrep}_run_@{i}.rds...\n")
			res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}/@{dataset}_subgroup_repeat_@{nrep}_run_@{i}.rds"))

			df_all = rbind(df_all, do.call(rbind, lapply(2:6, function(k) {
				df = get_stats(res, k = k)
				df = as.data.frame(df)
				df$method = rownames(df)
				df$i_simulate = i
				best_k = suggest_best_k(res)[rownames(df), ]
				df$best_k = best_k$best_k
				df$best_k_sig = best_k[, ncol(best_k)]
				df$running_time = sapply(res@list[rownames(df)], function(x) x@running_time)
				rownames(df) = NULL
				df
			})))
		}
		x = tapply(df_all[[stat]], paste0(df_all$method, df_all$k), mean)
		y = tapply(df_all[[stat]], paste0(df_all$method, df_all$k), sd)
		k = tapply(df_all$k, paste0(df_all$method, df_all$k), unique)
		method = tapply(df_all$method, paste0(df_all$method, df_all$k), unique)
		
		m = rbind(m, data.frame(x = x, y = y, k = k, nrep = nrep, method = method))
		df_all$nrep = nrep
		df_return = rbind(df_return, df_all)
	}
	m = as.data.frame(m)
	list(m = m, df_return = df_return)
}

# library(bsub)
# for(dataset in c("TCGA_GBM_by_row", "HSMM_single_cell_by_row", "TCGA_GBM_by_column", "HSMM_single_cell_by_column")) {
# 	for(stat in c("1-PAC", "mean_silhouette", "concordance")) {
# 		bsub_chunk(name = qq("nrep_@{dataset}_@{stat}"), hour = 5, memory = 5, 
# 			variables = c("dataset", "stat", "get_values"), packages = c("cola", "GetoptLong"),
# 		{
# 			m = get_values(stat)
# 			saveRDS(m, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}_subgroup_repeat_@{stat}.rds"))
# 		})
# 	}
# }

# dataset = "TCGA_GBM_by_row"
library(ggplot2)
library(scales)
library(ggrepel)
library(cowplot)
for(stat in c("1-PAC", "mean_silhouette", "concordance")) {
	file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}_subgroup_repeat_@{stat}.rds")
	lt = readRDS(file)
	m = as.data.frame(lt$m)
	m$x = as.numeric(m$x)
	m$y = as.numeric(m$y)
	m$cv = m$y/m$x
	df_all = lt$df_return

	pdf(qq("@{BASE_DIR}/image/@{dataset}_subgroup_repeat_@{stat}_mean_vs_sd.pdf"), width = 14, height = 4)
	p = ggplot(m, aes(x = x, y = y)) +
		geom_point(aes(col = factor(k))) + geom_smooth(se = FALSE, method = "loess", col = "grey", lwd = 1) +
		facet_wrap(~ factor(paste0(nrep, " samplings"), levels = c("25 samplings", "50 samplings", "100 samplings", "200 samplings")), nrow = 1) +
		xlab(qq("mean @{stat}")) + ylab("Standard deviation") + labs(col = "k")
	subset = m[order(m$y, decreasing = TRUE)[1:8], ]
	p = p + geom_text_repel(data = subset, aes(x = x, y = y, label = method))
	print(p)
	dev.off()

	if(stat == "1-PAC") {
		p = ggplot(m, aes(x = x, y = y)) +
			geom_point(aes(col = factor(k))) + geom_smooth(se = FALSE, method = "loess", col = "grey", lwd = 1) +
			facet_wrap(~ factor(paste0(nrep, " samplings"), levels = c("25 samplings", "50 samplings", "100 samplings", "200 samplings")), nrow = 2) +
			xlab(qq("mean @{stat}")) + ylab("Standard deviation") + labs(col = "k")
		subset = m[order(m$y, decreasing = TRUE)[1:8], ]
		p = p + geom_text_repel(data = subset, aes(x = x, y = y, label = method))
		print(p)

		if(dataset == "TCGA_GBM_by_row") {
			m = m[!(m$method == "ATC:mclust" & m$k == 2), ]
			p = ggplot(m, aes(x = x, y = y)) +
				geom_point(aes(col = factor(k))) + geom_smooth(se = FALSE, method = "loess", col = "grey", lwd = 1) +
				facet_wrap(~ factor(paste0(nrep, " samplings"), levels = c("25 samplings", "50 samplings", "100 samplings", "200 samplings")), nrow = 2) +
				xlab(qq("mean @{stat}")) + ylab("Standard deviation") + labs(col = "k")
		}
		figure_a = p
		# mfoo = data.frame(x = df_all[df_all$nrep == 25, stat], y = df_all[df_all$nrep == 200, stat], k = df_all[df_all$nrep == 25, "k"])
		# figure_b = ggplot(mfoo, aes(x = x, y = y)) + geom_point(aes(col = factor(k))) + xlab("25 samplings") + ylab("200 samplings") + labs(col = "k")
	}

	# pdf(qq("@{BASE_DIR}/image/@{dataset}_subgroup_repeat_@{stat}_nrep_compare.pdf"), width = 20, height = 6)
	# df = data.frame(x = df_all[df_all$nrep == 25, stat], y = df_all[df_all$nrep == 50, stat], k = df_all[df_all$nrep == 25, "k"])
	# p1 = ggplot(df, aes(x = x, y = y, col = factor(k))) + geom_point() + ggtitle(qq("@{stat}")) + xlab("25 samplings") + ylab("50 samplings") + labs(col = "k")
	
	# df = data.frame(x = df_all[df_all$nrep == 50, stat], y = df_all[df_all$nrep == 100, stat], k = df_all[df_all$nrep == 50, "k"])
	# p2 = ggplot(df, aes(x = x, y = y, col = factor(k))) + geom_point() + ggtitle(qq("@{stat}")) + xlab("50 samplings") + ylab("100 samplings") + labs(col = "k")
	
	# df = data.frame(x = df_all[df_all$nrep == 100, stat], y = df_all[df_all$nrep == 200, stat], k = df_all[df_all$nrep == 100, "k"])
	# p3 = ggplot(df, aes(x = x, y = y, col = factor(k))) + geom_point() + ggtitle(qq("@{stat}")) + xlab("100 samplings") + ylab("200 samplings") + labs(col = "k")
	
	# df = data.frame(x = df_all[df_all$nrep == 25, stat], y = df_all[df_all$nrep == 200, stat], k = df_all[df_all$nrep == 25, "k"])
	# p4 = ggplot(df, aes(x = x, y = y, col = factor(k))) + geom_point() + ggtitle(qq("@{stat}")) + xlab("25 samplings") + ylab("200 samplings") + labs(col = "k")
	
	# library(cowplot)
	# print(plot_grid(p1, p2, p3, p4, nrow = 1, align = "h"))
	# dev.off()
}


################
## running time
##############
ds = gsub('_by.*$', '', dataset)
file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{ds}_by_row_subgroup_repeat_1-PAC.rds")
lt = readRDS(file)
df_all_row = lt$df_return
df_all_row$by = "row"
file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{ds}_by_column_subgroup_repeat_1-PAC.rds")
lt = readRDS(file)
df_all_column = lt$df_return
df_all_column$by = "column"
df_all = rbind(df_all_row, df_all_column)
df_all$pm = gsub("^.*:", "", df_all$method)
df_all$pm = factor(df_all$pm, levels = c("hclust", "kmeans", "pam", "mclust", "skmeans", "NMF"))
df_all = df_all[df_all$k == 2, ]

library(ggplot2)
pdf(qq("@{BASE_DIR}/image/@{ds}_row_column_running_time.pdf"), width = 14, height = 8)
p1 = ggplot(df_all, aes(x = pm, y = running_time, col = by)) +
	geom_boxplot() +
	facet_wrap( ~ nrep, nrow = 1) +
	labs(x = "Partition methods", y = "Running time (sec)", col = "Sample by")

l = df_all$by == "row"
mean_running_time_row = tapply(df_all$running_time[l], paste0(df_all$pm[l], ":", df_all$nrep[l]), mean) 
l = df_all$by == "column"
mean_running_time_col = tapply(df_all$running_time[l], paste0(df_all$pm[l], ":", df_all$nrep[l]), mean) 
mean_running_time_col = mean_running_time_col[names(mean_running_time_row)]
ratio = mean_running_time_col/mean_running_time_row
df_mean_running_time = data.frame(
	pm = gsub(":.*$", "", names(ratio)),
	nrep = as.numeric(gsub("^.*:", "", names(ratio))),
	ratio = ratio
)
df_mean_running_time$pm = factor(df_mean_running_time$pm, levels = c("hclust", "kmeans", "pam", "mclust", "skmeans", "NMF"))
p2 = ggplot(df_mean_running_time, aes(x = pm, y = ratio)) +
	geom_bar(stat = "identity") +
	facet_wrap( ~ nrep, nrow = 1) +
	labs(x = "Partition methods", y = "Ratio of runnign time as column/row")

print(plot_grid(p1, p2, nrow = 2, labels = c("A", "B")))
dev.off()



foo = tapply(df_all$running_time, paste(df_all$nrep, df_all$pm, df_all$by), mean)
mean(foo[grepl("column", names(foo))]/foo[grepl("row", names(foo))])
mean(foo[grepl("NMF.*column", names(foo))]/foo[grepl("NMF.*row", names(foo))])

###########################################
## concordance 
###########################################

# we need to compare partitions between different nrep, thus, we need
# to adjust the class labels
res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}/@{dataset}_subgroup_repeat_25_run_1.rds"))
all_methods = names(res@list)

## read class id from all results
# initialize cl
cl = list()
for(method in all_methods) {
	cl[[method]] = list()
	for(k in 2:6) {
		cl[[method]][[as.character(k)]] = list()
		for(nrep in c(25, 50, 100, 200)) {
			cl[[method]][[as.character(k)]][[as.character(nrep)]] = list()
			for(i in 1:100) {
				cl[[method]][[as.character(k)]][[as.character(nrep)]][[i]] = NA
			}
		}
	}
}
# read original class id to cl
for(nrep in c(25, 50, 100, 200)) {
	for(i in 1:100) {
		qqcat("loading @{dataset}_subgroup_repeat_@{nrep}_run_@{i}.rds...\n")
		res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}/@{dataset}_subgroup_repeat_@{nrep}_run_@{i}.rds"))
		for(method in all_methods) {
			for(k in 2:6) {
				foo = get_classes(res[method], k = k)[, 1]
				attr(foo, "stat") = get_stats(res[method], k = k)
				cl[[method]][[as.character(k)]][[as.character(nrep)]][[i]] = foo
			}
		}
	}
}
saveRDS(cl, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}_subgroup_repeat_all_class_ids.rds"))

## for each method, each k, adjust class labels for different nrep

## consensus class id from 100 runs, weighted by mean silhouette (cross sample) in each run
library(clue)
make_global_class_ids = function(cl, method, k) {
	qqcat("formating membership matrix for @{method}, k = @{k}\n")
	
	cl_list = do.call("c", cl[[method]][[as.character(k)]])
	
	partition_list = lapply(cl_list, as.cl_hard_partition)
	partition_list = cl_ensemble(list = partition_list)
	partition_consensus = cl_consensus(partition_list)

	class_ids = as.vector(cl_class_ids(partition_consensus))
	
	return(class_ids)
}

global_cl = list()
for(method in all_methods) {
	global_cl[[method]] = list()
	for(k in 2:6) {
		global_cl[[method]][[as.character(k)]] = make_global_class_ids(cl, method, k)
	}
}
## adjust class IDs in `cl`, for each k
cl2 = cl
for(method in all_methods) {
	for(k in 2:6) {
		qqcat("relabel class IDs for @{method}, k = @{k}\n")
		for(nrep in c(25, 50, 100, 200)) {
			for(i in 1:100) {
				cl2[[method]][[as.character(k)]][[as.character(nrep)]][[i]] = relabel_class(cl[[method]][[as.character(k)]][[as.character(nrep)]][[i]], global_cl[[method]][[as.character(k)]], return_map = FALSE)
			}
		}
	}
}
saveRDS(cl2, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}_subgroup_repeat_all_class_ids2.rds"))

df = data.frame(method = character(0), k = numeric(0), nrep = numeric(0), p = numeric(0))
for(method in all_methods) {
	for(k in 2:6) {
		for(nrep in c(25, 50, 100, 200)) {

			lt1 = cl2[[method]][[as.character(k)]][[as.character(nrep)]]
			pp = numeric(length(lt1)*(length(lt1) - 1)/2)
			ii = 0		
			for(i in seq_along(lt1)[-1]) {
				for(j in 1:i) {
					c1 = lt1[[i]]
					c2 = lt1[[j]]
					ii = ii + 1
					pp[ ii ] = sum(c1 == c2)/length(c1)
				}
			}
			df = rbind(df, data.frame(method = method, k = k, nrep = nrep, p = mean(pp), stringsAsFactors = FALSE))
		}
	}
}
# for a method, a nrep, a k, the mean concordance between the 100 runs
pdf(qq("@{BASE_DIR}/image/@{dataset}_subgroup_repeat_concordance.pdf"), width = 14, height = 8)
f_scale = function(p) 1 - (1-p)^0.5
p = ggplot() + 
    geom_line(data = df, aes(x=log2(nrep/25)+1, y=f_scale(p), col = factor(k)), stat="identity") + 
    facet_grid(gsub(":.*$", "", method) ~ gsub("^.*:", "", method)) + 
    scale_y_continuous(breaks = f_scale(c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1)), 
    	labels = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1)) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(25, 50, 100, 200)) +
    ylab("corcodance") + labs(col = "k", x = "Number of samplings")
print(p)
dev.off()

figure_c = ggplot() + 
    geom_line(data = df[df$method == "MAD:hclust", ], aes(x=log2(nrep/25)+1, y=f_scale(p), col = factor(k)), stat="identity") + 
    scale_y_continuous(breaks = f_scale(c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1)), 
    	labels = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1)) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(25, 50, 100, 200)) +
    ylab("corcodance") + ggtitle("MAD:hclust") + labs(col = "k", x = "Number of samplings")


df2 = data.frame(method = character(0), k = numeric(0), p = numeric(0), pac_25 = numeric(0), pac_200 = numeric(0))
for(method in all_methods) {
	for(k in 2:6) {

		lt1 = cl2[[method]][[as.character(k)]][["25"]]
		lt2 = cl2[[method]][[as.character(k)]][["200"]]

		lt12 = cl[[method]][[as.character(k)]][["25"]]
		lt22 = cl[[method]][[as.character(k)]][["200"]]

		pp = numeric(length(lt1)*length(lt2))
		ii = 0	
		for(i in seq_along(lt1)) {
			for(j in seq_along(lt2)) {
				c1 = lt1[[i]]
				c2 = lt2[[j]]
				ii = ii + 1
				pp[ ii ] = sum(c1 == c2)/length(c1)
			}
		}
		pac_25 = NULL
		pac_200 = NULL
		for(i in seq_along(lt12)) {
			pac_25[i] = attr(lt12[[i]], "stat")[, "1-PAC"]
		}
		for(i in seq_along(lt22)) {
			pac_200[i] = attr(lt22[[i]], "stat")[, "1-PAC"]
		}
		df2 = rbind(df2, data.frame(method = method, k = k, p = mean(pp), pac_25 = mean(pac_25), pac_200 = mean(pac_200), stringsAsFactors = FALSE))
		
	}
}
# concordance of partitions of 25 vs 200 to mean 1-PAC with nrep = 25
f_scale = function(p) 1 - (1-p)^0.5
pdf(qq("@{BASE_DIR}/image/@{dataset}_subgroup_repeat_25_vs_200.pdf"), width = 14, height = 8)
p = ggplot(df2, aes(x=factor(k), y=f_scale(p))) + 
    geom_bar(stat="identity", position=position_dodge()) + facet_wrap( ~ method, nrow = 1) + 
    facet_grid(gsub(":.*$", "", method) ~ gsub("^.*:", "", method)) + 
    scale_y_continuous(breaks = f_scale(c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1)), 
    	labels = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1), limits = c(0, 1)) +
    ylab("concordance between 25 and 200 samplings") + labs("x" = "k")
print(p)
dev.off()

figure_d = ggplot(df2, aes(x = pac_25, y = p, col = factor(k))) + geom_point() + ylab("concordance") + labs("x" = "1-PAC, 25 samplings", col = "k")

library(cowplot)
theme_set(theme_grey())

pdf(qq("@{BASE_DIR}/image/@{dataset}_figure6.pdf"), width = 9, height = 6)
p = plot_grid(figure_c, figure_d, labels = c("B", "C"), align = 'v', nrow = 2)
p = plot_grid(figure_a, p, labels = c('A', ''), ncol = 2, rel_widths = c(2, 1.2))
print(p)
dev.off()

library(ggforce)
pdf(qq("@{BASE_DIR}/image/@{dataset}_subgroup_repeat_25_vs_200_pac_vs_corcondance.pdf"), width = 10, height = 5)
print(plot_grid(ggplot(df2, aes(x = pac_25, y = p, col = factor(k))) + geom_point() + ylab("concordance between 25 and 200 samplings") + 
	                labs("x" = "1-PAC, 25 samplings", col = "k"),
          ggplot(df2, aes(x = pac_200, y = p, col = factor(k))) + geom_point() + ylab("concordance between 25 and 200 samplings") + 
              labs("x" = "1-PAC, 200 samplings", col = "k"),
          align = "h"))
dev.off()


# for(dataset in c("TCGA_GBM_by_row", "HSMM_single_cell_by_row", "TCGA_GBM_by_column", "HSMM_single_cell_by_column")) {
	# cmd = qq("module load R/3.6.0; Rscript /icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_nrep/test_nrep_plot.R --dataset @{dataset}")
	# dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/")
	# cmd = qq("perl /desktop-home/guz/project/development/ngspipeline2/bsub_single_line.pl --hour 5 --memory 10 --name @{dataset}_subgroup_repeat_plot --dir @{dir} --command '@{cmd}' --enforce")
	# system(cmd)
# }



library(GetoptLong)
library(cola)
library(matrixStats)
library(ConsensusClusterPlus)
library(circlize)
library(ComplexHeatmap)

dataset = "TCGA_GBM"
GetoptLong("dataset=s", "dataset")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"
root = "/desktop-home/guz"

rl = readRDS(qq("@{root}/project/development/cola_examples/@{dataset}/@{dataset}_subgroup.rds"))
m = get_matrix(rl)

test_cc = function(m2, pItem = 0.8, pFeature = 1, clusterAlg = "km", n = 100) {
	lt = list()
	for(i in 1:n) {
		qqcat("running CC with @{clusterAlg}, with pItem @{pItem}, pFeature @{pFeature}, @{i}/@{n}\n")
		pdf(NULL)
		results = ConsensusClusterPlus(m2,maxK=6,reps=100,pItem=pItem,pFeature=pFeature,clusterAlg=clusterAlg)
		dev.off()
		
		cm_list = list()
		cl_list = list()
		for(k in 2:6) {
			cm_list[[k - 1]] = results[[k]]$consensusMatrix
			cl_list[[k - 1]] = results[[k]]$consensusClass
		}

		lt[[i]] = list(consensusMatrix = cm_list,
			           consensusClass = cl_list)
	}
	return(lt)
}

v1 = rowSds(m); od1 = order(v1, decreasing = TRUE)[1:2000]
m2 = t(scale(t(m[od1, ])))

set.seed(12345)
for(clusterAlg in c("km", "pam", "hc")) {
	lt = test_cc(m2, pItem = 0.8, pFeature = 1, clusterAlg = clusterAlg)
	saveRDS(lt, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}_row_or_column_resampling_top2k_sd_@{clusterAlg}_pItem0.8_pFeature1.rds"))
	lt = test_cc(m2, pItem = 1, pFeature = 0.8, clusterAlg = clusterAlg)
	saveRDS(lt, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}_row_or_column_resampling_top2k_sd_@{clusterAlg}_pItem1_pFeature0.8.rds"))
	lt = test_cc(m2, pItem = 0.8, pFeature = 0.8, clusterAlg = clusterAlg)
	saveRDS(lt, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}_row_or_column_resampling_top2k_sd_@{clusterAlg}_pItem0.8_pFeature0.8.rds"))
}

v2 = ATC(m); od2 = order(v2, decreasing = TRUE)[1:2000]
m2 = t(scale(t(m[od2, ])))
set.seed(12345)
for(clusterAlg in c("km", "pam", "hc")) {
	lt = test_cc(m2, pItem = 0.8, pFeature = 1, clusterAlg = clusterAlg)
	saveRDS(lt, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}_row_or_column_resampling_top2k_ATC_@{clusterAlg}_pItem0.8_pFeature1.rds"))
	lt = test_cc(m2, pItem = 1, pFeature = 0.8, clusterAlg = clusterAlg)
	saveRDS(lt, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}_row_or_column_resampling_top2k_ATC_@{clusterAlg}_pItem1_pFeature0.8.rds"))
	lt = test_cc(m2, pItem = 0.8, pFeature = 0.8, clusterAlg = clusterAlg)
	saveRDS(lt, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}_row_or_column_resampling_top2k_ATC_@{clusterAlg}_pItem0.8_pFeature0.8.rds"))
}

# for(dataset in c("Golub_leukemia", "HSMM_single_cell", "MCF10CA_scRNAseq", "Ritz_ALL", "TCGA_GBM")) {
# 	cmd = qq("module load R/3.3.1; Rscript /icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/row_or_column_resampling.R --dataset @{dataset}")
# 	cmd = qq("perl /desktop-home/guz/project/development/ngspipeline2/bsub_single_line.pl --hour 50 --memory 20 --name @{dataset}_row_or_column_resampling --command '@{cmd}' --enforce")
# 	system(cmd)
# }

# dataset in c("Golub_leukemia", "HSMM_single_cell", "MCF10CA_scRNAseq", "Ritz_ALL", "TCGA_GBM")

# dataset = "Golub_leukemia"
get_values = function(top_value_method, clusterAlg, pItem = 1, pFeature = 1) {
	lt = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}_row_or_column_resampling_top2k_@{top_value_method}_@{clusterAlg}_pItem@{pItem}_pFeature@{pFeature}.rds"))

	pac = NULL
	pac_overall = NULL
	for(k in 2:6) {
		pac = cbind(pac, sapply(lt, function(x) PAC(x$consensusMatrix[[k - 1]])))
		partition_list = lapply(lt, function(x) x$consensusClass[[k - 1]])
		pac_overall[ k - 1 ] = cal_consensus(partition_list, k)$stat["PAC"]
	}
	colnames(pac) = paste0("k=", 2:6)
	names(pac_overall) = paste0("k=", 2:6)

	return(list(pac = pac, 
		pac_overall = pac_overall, 
		best_k = apply(pac, 1, function(x) {
			ind = which.min(x)
			ind + 1
		}),
		best_k_pac = apply(pac, 1, min)))
}


lt = list(
	sd_hc_pItem0.8_pFeature1 = get_values("sd", "hc", 0.8, 1),
	sd_hc_pItem1_pFeature0.8 = get_values("sd", "hc", 1, 0.8),
	sd_km_pItem0.8_pFeature1 = get_values("sd", "km", 0.8, 1),
	sd_km_pItem1_pFeature0.8 = get_values("sd", "km", 1, 0.8),
	sd_pam_pItem0.8_pFeature1 = get_values("sd", "pam", 0.8, 1),
	sd_pam_pItem1_pFeature0.8 = get_values("sd", "pam", 1, 0.8),

	ATC_hc_pItem0.8_pFeature1 = get_values("ATC", "hc", 0.8, 1),
	ATC_hc_pItem1_pFeature0.8 = get_values("ATC", "hc", 1, 0.8),
	ATC_km_pItem0.8_pFeature1 = get_values("ATC", "km", 0.8, 1),
	ATC_km_pItem1_pFeature0.8 = get_values("ATC", "km", 1, 0.8),
	ATC_pam_pItem0.8_pFeature1 = get_values("ATC", "pam", 0.8, 1),
	ATC_pam_pItem1_pFeature0.8 = get_values("ATC", "pam", 1, 0.8)
)

df_all = data.frame(sample_by = character(0),
	method = character(0),
	k = numeric(0),
	pac = numeric(0))
for(i in seq_along(lt)) {
	for(k in 1:5) {
		df = data.frame(sample_by = ifelse(i %% 2 == 1, "by sample", "by gene"),
						method = gsub("_pItem.*$", "", names(lt)[i]),
			            k = k + 1,
			            pac = lt[[i]]$pac[, k],
			            stringsAsFactors = FALSE)
		df_all = rbind(df, df_all)
	}
}

df_all$sample_by = factor(df_all$sample_by, levels = c("by gene", "by sample"))
df_all$method = factor(df_all$method, levels = c("sd_hc", "sd_km", "sd_pam", "ATC_hc", "ATC_km", "ATC_pam"))
pdf(qq("@{BASE_DIR}/image/@{dataset}_row_or_column_resampling_boxplot.pdf"), width = 14, height = 5)
library(ggplot2)
library(scales)
ggplot(df_all, aes(x=factor(k), y=pac, fill=factor(sample_by))) + 
    geom_boxplot() + facet_wrap( ~ method, nrow = 1)
dev.off()


### consistency of best_k significant
df_all = data.frame(sample_by = character(0),
	method = character(0),
	k = numeric(0),
	p = numeric(0))
for(i in seq_along(lt)) {
	tb = table(lt[[i]]$best_k[lt[[i]]$best_k_pac <= 0.15])
	v = numeric(5)
	names(v) = 2:6
	v[names(tb)] = tb
	df = data.frame(sample_by = ifelse(i %% 2 == 1, "by sample", "by gene"),
					method = gsub("_pItem.*$", "", names(lt)[i]),
		            k = 2:6,
		            p = v/100,
		            stringsAsFactors = FALSE)
	df_all = rbind(df, df_all)
}

df_all$sample_by = factor(df_all$sample_by, levels = c("by gene", "by sample"))
df_all$method = factor(df_all$method, levels = c("sd_hc", "sd_km", "sd_pam", "ATC_hc", "ATC_km", "ATC_pam"))


pdf(qq("@{BASE_DIR}/image/@{dataset}_row_or_column_resampling_best_k_pac_0.15.pdf"), width = 14, height = 5)
ggplot(df_all, aes(x=method, y=p, fill=factor(k))) + 
    geom_bar( stat="identity") + 
    facet_wrap(~ sample_by)
dev.off()


### consistency of partition
df_all = data.frame(sample_by = character(0),
	method = character(0),
	k = numeric(0),
	pac_overall = numeric(0))
for(i in seq_along(lt)) {
	v = lt[[i]]$pac_overall
	df = data.frame(sample_by = ifelse(i %% 2 == 1, "by sample", "by gene"),
					method = gsub("_pItem.*$", "", names(lt)[i]),
		            k = 2:6,
		            pac_overall = v,
		            stringsAsFactors = FALSE)
	df_all = rbind(df, df_all)
}

df_all$pac_overall[df_all$pac_overall == 0] = 0.000001
df_all$sample_by = factor(df_all$sample_by, levels = c("by gene", "by sample"))
df_all$method = factor(df_all$method, levels = c("sd_hc", "sd_km", "sd_pam", "ATC_hc", "ATC_km", "ATC_pam"))

pdf(qq("@{BASE_DIR}/image/@{dataset}_row_or_column_resampling_pac_overall.pdf"), width = 14, height = 5)
ggplot(df_all, aes(x=k, y=pac_overall, fill=factor(sample_by))) + 
    geom_bar(stat="identity", position="dodge") + 
    scale_y_continuous(trans = sqrt_trans(),
    	breaks = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)) +
    facet_wrap( ~ method, nrow = 1)
dev.off()


get_membership_mat = function(top_value_method, clusterAlg, pItem = 1, pFeature = 1) {
	lt = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}_row_or_column_resampling_top2k_@{top_value_method}_@{clusterAlg}_pItem@{pItem}_pFeature@{pFeature}.rds"))

	ml = list()
	for(k in 2:6) {
		ml[[k - 1]] = lapply(lt, function(x) x$consensusClass[[k - 1]])
	}
	ml
}


lt = list(
	sd_hc_pItem0.8_pFeature1 = get_membership_mat("sd", "hc", 0.8, 1),
	sd_hc_pItem1_pFeature0.8 = get_membership_mat("sd", "hc", 1, 0.8),
	sd_km_pItem0.8_pFeature1 = get_membership_mat("sd", "km", 0.8, 1),
	sd_km_pItem1_pFeature0.8 = get_membership_mat("sd", "km", 1, 0.8),
	sd_pam_pItem0.8_pFeature1 = get_membership_mat("sd", "pam", 0.8, 1),
	sd_pam_pItem1_pFeature0.8 = get_membership_mat("sd", "pam", 1, 0.8),

	ATC_hc_pItem0.8_pFeature1 = get_membership_mat("ATC", "hc", 0.8, 1),
	ATC_hc_pItem1_pFeature0.8 = get_membership_mat("ATC", "hc", 1, 0.8),
	ATC_km_pItem0.8_pFeature1 = get_membership_mat("ATC", "km", 0.8, 1),
	ATC_km_pItem1_pFeature0.8 = get_membership_mat("ATC", "km", 1, 0.8),
	ATC_pam_pItem0.8_pFeature1 = get_membership_mat("ATC", "pam", 0.8, 1),
	ATC_pam_pItem1_pFeature0.8 = get_membership_mat("ATC", "pam", 1, 0.8)
)

make_membership_mat = function(lt, k) {
	lt2 = lapply(lt, function(x) x[[k - 1]])
	lt3 = do.call(c, lt2)
	nm = rep(names(lt2), each = length(lt2[[1]]))
	res = cal_consensus(lt3, k)
	m = res$membership_each

	lapply(split(seq_along(nm), nm), function(x) t(m[, x]))
}


pdf(qq("@{BASE_DIR}/image/@{dataset}_row_or_column_resampling_membership.pdf"), width = 5, height = 10)
for(k in 2:6) {
	ml = make_membership_mat(lt, k = k)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nr = 6, nc = 2)))
	for(i in seq_along(ml)) {
		pushViewport(viewport(layout.pos.row = ceiling(i / 2), layout.pos.col = (i + 1) %% 2 + 1))
		if((i + 1) %% 2 + 1 == 1) {
			od = hclust(dist(t(ml[[i]])))$order
		}
		draw(Heatmap(ml[[i]], col = cola:::brewer_pal_set2_col, use_raster = TRUE,
			show_heatmap_legend = FALSE, show_row_dend = FALSE,
			cluster_columns = FALSE, column_order = od), 
			newpage = FALSE)
		popViewport()
	}
	popViewport()
}
dev.off()


make_class_id = function(lt, k) {
	lt2 = lapply(lt, function(x) x[[k - 1]])
	lt3 = do.call(c, lt2)
	nm = rep(names(lt2), each = length(lt2[[1]]))
	res = cal_consensus(lt3, k)
	cl_global = res$class_ids

	do.call("rbind", lapply(lt2, function(x) {
		cl = cal_consensus(x, k)$class_ids
		# adjust according to cl_global
		map = relabel_class(cl, cl_global, full_set = 1:k)
		as.integer(as.numeric(map[as.character(cl)]))
	}))
}

pdf(qq("@{BASE_DIR}/image/@{dataset}_row_or_column_resampling_class_id.pdf"), width = 5, height = 6)
for(k in 2:6) {
	m = make_class_id(lt, k = k)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nr = 6, nc = 1)))
	for(i in 1:6) {
		pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
		draw(Heatmap(m[2*i - c(1, 0), ], col = cola:::brewer_pal_set2_col,
			row_title = gsub("_pItem.*$", "", rownames(m)[2*i]),
			row_title_rot = 90,
			row_labels = c("by genes", "by samples"),
			cluster_rows = FALSE,
			show_heatmap_legend = FALSE, show_row_dend = FALSE, show_column_dend = FALSE), 
			newpage = FALSE)
		popViewport()
	}
	popViewport()
}
dev.off()


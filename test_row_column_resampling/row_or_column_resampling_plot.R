library(GetoptLong)
library(cola)
library(matrixStats)
library(ConsensusClusterPlus)
library(circlize)
library(ComplexHeatmap)

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"

#dataset in c("Golub_leukemia", "HSMM_single_cell", "MCF10CA_scRNAseq", "Ritz_ALL", "TCGA_GBM"

dataset = "Golub_leukemia"
GetoptLong("dataset=s", "dataset")

df_all = data.frame(k = numeric(0), `1-PAC` = numeric(0), method = character(0), sample_by = character(0),
	check.names = FALSE)
for(i in 1:100) {
	qqcat("loading @{dataset}, @{i}/100, ...\n")
	res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}/@{dataset}_resampling_by_row_@{i}.rds"))
	for(k in 2:6) {	
		stat_df = as.data.frame(get_stats(res, k = k))
		stat_df$method = rownames(stat_df)
		stat_df$sample_by = "by row"
		df_all = rbind(df_all, stat_df[, c("k", "1-PAC", "method", "sample_by")])
	}
	res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}/@{dataset}_resampling_by_column_@{i}.rds"))
	for(k in 2:6) {	
		stat_df = as.data.frame(get_stats(res, k = k))
		stat_df$method = rownames(stat_df)
		stat_df$sample_by = "by column"
		df_all = rbind(df_all, stat_df[, c("k", "1-PAC", "method", "sample_by")])
	}
}

df_all$sample_by = factor(df_all$sample_by, levels = c("by row", "by column"))
df_all$method = factor(df_all$method, levels = c("SD:hclust", "SD:skmeans", "ATC:hclust", "ATC:skmeans"))


pdf(qq("@{BASE_DIR}/image/@{dataset}_row_or_column_resampling_boxplot.pdf"), width = 14, height = 4)
library(ggplot2)
library(scales)
p = ggplot(df_all, aes(x=factor(k), y=`1-PAC`, fill=sample_by)) + 
    geom_boxplot() + facet_wrap( ~ method, nrow = 1) + ylim(0, 1) +
    labs(x = "k", fill = "Sample by")
print(p)
dev.off()

figure_a = p

# median difference
df_all1 = df_all[df_all$sample_by == "by row", ]
x1 = tapply(df_all1$`1-PAC`, paste(df_all1$method, df_all1$k, sep = "-"), median)
df_all2 = df_all[df_all$sample_by == "by column", ]
x2 = tapply(df_all2$`1-PAC`, paste(df_all2$method, df_all2$k, sep = "-"), median)
diff = x1 - x2

diff_df = data.frame(method = gsub("-\\d$", "", names(diff)),
	k = as.numeric(gsub("^.*-", "", names(diff))),
	diff = diff)
diff_df$method = factor(diff_df$method, levels = c("SD:hclust", "SD:skmeans", "ATC:hclust", "ATC:skmeans"))

pdf(qq("@{BASE_DIR}/image/@{dataset}_row_or_column_resampling_median_diff.pdf"), width = 6, height = 4)
library(ggplot2)
library(scales)
p = ggplot(diff_df, aes(x=factor(k), y=diff)) + 
    geom_bar(stat="identity") + facet_wrap( ~ method, nrow = 1) + 
    ylab("1-PAC difference") +
    ggtitle("1-PAC difference (by row - by column)") +
    labs(x = "k")
print(p)
dev.off()

figure_b = p

#######################################
## global class id for each k
library(clue)
lt = lapply(2:6, function(x) lapply(1:4, function(y) list()))
lt_col = lt_row = lt
for(i in 1:100) {
	qqcat("loading @{dataset}, @{i}/100, ...\n")
	res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}/@{dataset}_resampling_by_row_@{i}.rds"))
	for(k in 2:6) {	
		for(j in seq_along(res@list)) {
			lt[[k - 1]][[j]][[i*2-1]] = get_classes(res@list[[j]], k = k)[, 1]
			lt_row[[k - 1]][[j]][[i]] = get_classes(res@list[[j]], k = k)[, 1]
		}
		names(lt[[k-1]]) = names(res@list)
		names(lt_row[[k-1]]) = names(res@list)
	}

	res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}/@{dataset}_resampling_by_column_@{i}.rds"))
	for(k in 2:6) {	
		for(j in seq_along(res@list)) {
			lt[[k - 1]][[j]][[i*2]] = get_classes(res@list[[j]], k = k)[, 1]
			lt_col[[k - 1]][[j]][[i]] = get_classes(res@list[[j]], k = k)[, 1]
		}
		names(lt[[k-1]]) = names(res@list)
		names(lt_col[[k-1]]) = names(res@list)
	}
}


# adjust to previous k
global_class_id = lapply(lt, function(x) {
	lapply(x, function(y) {
		y = lapply(y, as.cl_hard_partition)
		partition_list = cl_ensemble(list = y)
		partition_consensus = cl_consensus(partition_list)
		as.vector(cl_class_ids(partition_consensus))
	})
})

for(k in 3:6) {
	for(nm in names(global_class_id[[k-1]])) {
		global_class_id[[k - 1]][[nm]] = as.numeric(relabel_class(global_class_id[[k - 1]][[nm]], global_class_id[[k - 2]][[nm]], full_set = 1:k, return_map = FALSE))
	}
}


#######################################################
## how consistent the row sampling and column sampling
#######################################################
pdf(qq("@{BASE_DIR}/image/@{dataset}_row_or_column_resampling_heatmap.pdf"), width = 20, height = 16)
i = 0; j = 0
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 5)))
for(nm in c("SD:hclust", "SD:skmeans", "ATC:hclust", "ATC:skmeans")) {

	i = i + 1
	j = 0
	for(k in 2:6) {
		j = j + 1

		lt1 = lt_row[[k - 1]][[nm]]
		lt2 = lt_col[[k - 1]][[nm]]

		lt1 = lapply(lt1, function(x) relabel_class(x, global_class_id[[k - 1]][[nm]], return_map = FALSE))
		lt2 = lapply(lt2, function(x) relabel_class(x, global_class_id[[k - 1]][[nm]], return_map = FALSE))
		m1 = do.call("rbind", lt1)
		m2 = do.call("rbind", lt2)
		m = rbind(m1, m2)
		by = c(rep("by row", nrow(m1)), rep("by column", nrow(m2)))
		ht = Heatmap(m, 
			row_split = factor(by, levels = c("by row", "by column")), 
			column_split = factor(global_class_id[[k - 1]][[nm]], levels = c(1:k)),
			column_title = qq("@{nm}, @{k} subgroups"), show_row_dend = FALSE, show_column_dend = FALSE,
			col = cola:::brewer_pal_set2_col, name = "subgroup", use_raster = TRUE,
			cluster_column_slices = FALSE, cluster_row_slices = FALSE)
		pushViewport(viewport(layout.pos.row = i, layout.pos.col = j))
		draw(ht, newpage = FALSE)
		popViewport()
	}
}
popViewport()
dev.off()

# or compare any pair of partitions
# 1. only for by row 
# 2. only for by column
# 3. by row ~ by column
library(parallel)
p_between = p_row = p_col = NULL
for(nm in c("SD:hclust", "SD:skmeans", "ATC:hclust", "ATC:skmeans")) {

	pl = mclapply(2:6, function(k) {

		lt1 = lt_row[[k - 1]][[nm]]
		lt2 = lt_col[[k - 1]][[nm]]

		lt1 = lapply(lt1, function(x) relabel_class(x, global_class_id[[k - 1]][[nm]], return_map = FALSE))
		lt2 = lapply(lt2, function(x) relabel_class(x, global_class_id[[k - 1]][[nm]], return_map = FALSE))

		qqcat("calculte mean concordance for @{nm}_k@{k}\n")

		pp = numeric(length(lt1)*length(lt2))	
		ii = 0	
		for(i in seq_along(lt1)) {
			for(j in seq_along(lt2)) {
				cl1 = lt1[[i]]
				cl2 = lt2[[j]]
				ii = ii + 1
				# pp[ ii ] = sum(cl1 == cl2)/length(cl1)
				pp[ ii ] = mclust::adjustedRandIndex(cl1, cl2)
			}
		}
		p_between = structure(names = qq("@{nm}_k@{k}"), mean(pp))
		
		pp = numeric(length(lt1)*(length(lt1) - 1)/2)
		ii = 0		
		for(i in seq_along(lt1)[-1]) {
			for(j in 1:i) {
				cl1 = lt1[[i]]
				cl2 = lt1[[j]]
				ii = ii + 1
				# pp[ ii ] = sum(cl1 == cl2)/length(cl1)
				pp[ ii ] = mclust::adjustedRandIndex(cl1, cl2)
			}
		}
		p_row = structure(names = qq("@{nm}_k@{k}"), mean(pp))

		pp = numeric(length(lt2)*(length(lt2) - 1)/2)
		ii = 0		
		for(i in seq_along(lt2)[-1]) {
			for(j in 1:i) {
				cl1 = lt2[[i]]
				cl2 = lt2[[j]]
				ii = ii + 1
				# pp[ ii ] = sum(cl1 == cl2)/length(cl1)
				pp[ ii ] = mclust::adjustedRandIndex(cl1, cl2)
			}
		}
		p_col = structure(names = qq("@{nm}_k@{k}"), mean(pp))
		
		c(p_between, p_row, p_col)
	}, mc.cores = 5)
	p_between = c(p_between, sapply(pl, "[", 1))
	p_row = c(p_row, sapply(pl, "[", 2))
	p_col = c(p_col, sapply(pl, "[", 3))
}


df1 = data.frame(
	method = gsub("_k\\d$", "", names(p_between)),
	k = gsub("^.*_k", "", names(p_between)),
	p = p_between,
	compare = "row ~ column",
	stringsAsFactors = FALSE)
df2 = data.frame(
	method = gsub("_k\\d$", "", names(p_col)),
	k = gsub("^.*_k", "", names(p_col)),
	p = p_col,
	compare = "column",
	stringsAsFactors = FALSE)
df3 = data.frame(
	method = gsub("_k\\d$", "", names(p_row)),
	k = gsub("^.*_k", "", names(p_row)),
	p = p_row,
	compare = "row",
	stringsAsFactors = FALSE)
df = rbind(df1, df2, df3)
df$method = factor(df$method, levels = c("SD:hclust", "SD:skmeans", "ATC:hclust", "ATC:skmeans"))
df$compare = factor(df$compare, levels = c("row ~ column", "column", "row"))

f_scale = function(p) 1 - (1-p)^0.5
library(ggplot2)
pdf(qq("@{BASE_DIR}/image/@{dataset}_row_or_column_resampling_corcodance.pdf"), width = 14, height = 4)
p = ggplot(df, aes(x=factor(k), y=f_scale(p), fill = compare)) + 
    geom_bar(stat="identity", position=position_dodge()) + facet_wrap( ~ method, nrow = 1) + 
    geom_hline(yintercept = f_scale(0.9), color = "grey", lty = 2) +
    scale_y_continuous(breaks = f_scale(c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1)), 
    	labels = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1), limits = c(0, 1)) +
    ylab("Corcodance between row/column sampling") +
    labs(x = "k", fill = "Comparison")
print(p)
dev.off()

p = ggplot(df[df$compare == "row ~ column", ], aes(x=factor(k), y=p)) + 
    geom_bar(stat="identity", position=position_dodge()) + facet_wrap( ~ method, nrow = 1) + 
    scale_y_continuous(limits = c(0, 1)) +
    ylab("Corcodance") +
    ggtitle("Corcodance between row/column sampling") +
    labs(x = "k")

figure_c = p

library(cowplot)
theme_set(theme_grey())

pdf(qq("@{BASE_DIR}/image/@{dataset}_figure5.pdf"), width = 12, height = 6)
p = plot_grid(figure_a, 
		plot_grid(figure_b, figure_c, labels = c('B', 'C'), nrow = 1, align = "h"),
		align = "v", ncol = 1, labels = c("A", ""))
print(p)
dev.off()


# library(bsub)
# for(dataset in c("Golub_leukemia", "HSMM_single_cell", "MCF10CA_scRNAseq", "Ritz_ALL", "TCGA_GBM")) {
# 	bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_row_column_resampling/row_or_column_resampling_plot.R",
# 		name = qq("@{dataset}_row_or_column_resampling_plot"), argv = qq("--dataset @{dataset}"),
# 		hour = 1, memory = 1)
# }



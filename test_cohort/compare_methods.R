

library(cola)
library(GetoptLong)
library(matrixStats)
BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"

####################################################
# project = "GDS"

GetoptLong("project=s", "project")

collect_all_stats = function(project) {
	if(project == "recount2") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
	} else if(project == "GDS") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
	}

	df = file.info(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id_list}", collapse = FALSE))
	id_list = basename(rownames(df)[df$isdir])

	lt = list()
	i = 0
	for(id in id_list) {
		i = i + 1
		qqcat("loading @{id} (@{i}/@{length(id_list)})...\n")
		oe = try({
			res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_cola_all.rds"))

			df = do.call(rbind, lapply(2:6, function(k) {
				df = get_stats(res, k = k)
				df = as.data.frame(df)
				df$method = rownames(df)
				df$id = id
				best_k = suggest_best_k(res)[rownames(df), ]
				df$best_k = best_k$best_k
				df$n_sample = ncol(res)
				
				if(is.null(get_anno(res))) {
					df$p_anno = NA
				} else {
					pmat = test_to_known_factors(res, k = k)
					pmat = pmat[rownames(df), , drop = FALSE]
					pmat = pmat[, -1, drop = FALSE]
					pmat = pmat[, -ncol(pmat), drop = FALSE]
					pmat = as.matrix(pmat)
					df$p_anno = rowMins(pmat)
				}
				rownames(df) = NULL
				df
			}))
		})

		if(!inherits(oe, "try-error")) {
			lt = c(lt, list(df))
		}
	}
	return(do.call("rbind", lt))
}

file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{project}_df_all.rds")
if(file.exists(file)) {
	df_all = readRDS(file)
} else {
	df_all = collect_all_stats(project)
	saveRDS(df_all, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{project}_df_all.rds"))
}

# density heatmap for the statistics
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(matrixStats)
make_plot = function(df_all, field) {
	lt = split(df_all, df_all$method)
	lt = lapply(lt, function(df) df[df$k == df$best_k, ])
	lt = lapply(lt, function(x) {
		x = x[, field]
		x[!is.na(x)]
	})

	anno = data.frame(
		good = ifelse(grepl("ATC|skmeans", names(lt)), "yes", "no"),
		bad = ifelse(grepl("hclust", names(lt)), "yes", "no"))
	anno_col = list(good = c("yes" = "red", "no" = "#EEEEEE"),
		            bad = c("yes" = "blue", "no" = "#EEEEEE"))

	ht = densityHeatmap(lt, column_order = order(sapply(lt, median)), ylab = field, ylim = c(min(unlist(lt), na.rm = TRUE), 1),
		bottom_annotation = HeatmapAnnotation(df = anno, col = anno_col, show_legend = FALSE, simple_anno_size = unit(2, "mm"), show_annotation_name = FALSE, gp = gpar(col = "white", lwd = 0.5)),
		column_names_rot = 45, column_title = qq("Density heatmap of @{field}"),
		column_names_gp = gpar(fontsize = 8))
	draw(ht)
	return(names(lt)[order(sapply(lt, median))])
}

image1 = grid.grabExpr({col_od <- make_plot(df_all, "1-PAC")})
image2 = grid.grabExpr(make_plot(df_all, "mean_silhouette"))
image3 = grid.grabExpr(make_plot(df_all, "concordance"))
pdf(qq("@{BASE_DIR}/image/@{project}_stat_best_k_density_distribution.pdf"), width = 20, height = 8)
grid.newpage()
pushViewport(viewport(x = 0, width = 1/3, just = "left"))
grid.draw(image1)
popViewport()
pushViewport(viewport(x = 1/3, width = 1/3, just = "left"))
grid.draw(image2)
popViewport()
pushViewport(viewport(x = 2/3, width = 1/3, just = "left"))
grid.draw(image3)
popViewport()
dev.off()

figure_b = image1

lt = split(df_all, df_all$method)
lt = lapply(lt, function(df) df[df$k == df$best_k, ])


# for each method, proportion of stable partitions for the best k
pct1 = sapply(sapply(lt, function(df) {
	structure(names = rownames(df), df$`1-PAC` >= 0.9)
}), function(x) sum(x, na.rm = TRUE)/length(x))

# proportion of best k
pct2 = sapply(sapply(lt, function(df) {
	best_k = df$best_k
	best_k[df$`1-PAC` < 0.9] = NA
	structure(names = rownames(df), best_k)
}), function(x) {
	tb = table(x)
	y = rep(0, 5)
	names(y) = 2:6
	y[names(tb)] = tb/length(x)
	y
})

library(reshape2)
library(ggplot2)

p2_df = melt(pct2)
colnames(p2_df) = c("k", "method", "p")

pdf(qq("@{BASE_DIR}/image/@{project}_best_k_prop.pdf"), width = 8, height = 8)

# ggplot(data.frame(p = pct1), aes(x = factor(names(pct1), levels = names(p1)[order(pct1)]), y = p)) + 
# 	geom_bar(stat = "identity") +
# 	theme(axis.text.x = element_text(angle = 45, hjust = 1))
	
pc = ggplot(p2_df, aes(x = factor(method, levels = col_od), y = p, fill = factor(k, levels = 6:2))) + 
	geom_bar(stat = "identity") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Partition methods", y = "Proportion", fill = "k") +
	ggtitle("Proportion of the best k")
print(pc)
dev.off()

figure_c = pc


library(cowplot)
theme_set(theme_grey())
pdf(qq("@{BASE_DIR}/image/@{project}_stat_ggparis.pdf"), width = 12, height = 4)
for(method in unique(df_all$method)) {
	l = df_all$method == method
	df2 = df_all[l, ]
	p1 = ggplot(df2, aes(x = mean_silhouette, y = concordance, color = paste0("k=", k))) +
		geom_point() + ggtitle(method) + labs(color = "Number of\npartitions")
	p2 = ggplot(df2, aes(x = mean_silhouette, y = `1-PAC`, color = paste0("k=", k))) +
		geom_point() + ggtitle(method) + labs(color = "Number of\npartitions")
	p3 = ggplot(df2, aes(x = concordance, y = `1-PAC`, color = paste0("k=", k))) +
		geom_point() + ggtitle(method) + labs(color = "Number of\npartitions")
	p = plot_grid(p1, p2, p3, nrow = 1)
	print(p)
}
dev.off()

foo = tapply(seq_len(nrow(df_all)), df_all[, c("k", "method")], function(ind) {
	cm = cor(df_all[ind, c("1-PAC", "mean_silhouette", "concordance")], method = "spearman")
	c("mean_silhouette\nvs\n1-PAC" = cm["mean_silhouette", "1-PAC"],
	  "concordance\nvs\n1-PAC" = cm["concordance", "1-PAC"],
	  "concordance\nvs\nmean_silhouette" = cm["concordance", "mean_silhouette"])
})

arr = array(0, dim = c(nrow(foo), ncol(foo), length(foo[1, 1][[1]])),
	dimnames = list(rownames(foo), colnames(foo), names(foo[1, 1][[1]])))
for(i in seq_len(nrow(foo))) {
	for(j in seq_len(ncol(foo))) {
		arr[i, j, ] = foo[i, j][[1]]
	}
}

library(reshape2)
df = melt(arr)
colnames(df) = c("k", "method", "comparison", "correlation")
pdf(qq("@{BASE_DIR}/image/@{project}_stat_correlation.pdf"), width = 7, height = 7)
p = ggplot(df, aes(x = comparison, y = correlation, fill = factor(k))) + 
	geom_boxplot() + ylim(c(0, 1)) + ggtitle("Correlation between statistics") +
	labs("x" = "Comparison", "y" = "Correlation", fill = "k")
print(p)
dev.off()

figure_a = p


####### dissimilarity between methods
library(clue)
collect_hclust = function(project, k) {
	if(project == "recount2") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
	} else if(project == "GDS") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
	}

	df = file.info(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id_list}", collapse = FALSE))
	id_list = basename(rownames(df)[df$isdir])

	lt = list()
	i = 0
	for(id in id_list) {
		i = i + 1
		qqcat("loading @{id} (@{i}/@{length(id_list)})...\n")
			
		res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_cola_all.rds"))

		pl = lapply(res@list, function(x) as.cl_partition(get_membership(x, k = k)))
		clen = cl_ensemble(list = pl)
		m_diss = cl_dissimilarity(clen, method = "euclidean")
		lt[[i]] = hclust(m_diss)
	}
	names(lt) = id_list
	return(lt)
}

pdf(qq("@{BASE_DIR}/image/@{project}_partition_similarity.pdf"), width = 7, height = 7)
k = 2
hl = collect_hclust(project, k = k)
hens = cl_ensemble(list = hl)
dend = cl_consensus(hens)
dend = as.dendrogram(hclust(dend$.Data))
library(ggdendro)

p = ggdendrogram(dend, rotate = TRUE) + scale_y_continuous(expand = c(0.01, 0.01), labels = NULL) +
    ggtitle(qq("Similarity of Partitions, k = @{k}"))
print(p)

label_color = c(`ATC:hclust` = "#66C2A5", `ATC:kmeans` = "#66C2A5", `ATC:pam` = "#66C2A5", 
`ATC:mclust` = "#66C2A5", `ATC:NMF` = "#66C2A5", `ATC:skmeans` = "#66C2A5", 
`CV:hclust` = "#FC8D62", `SD:hclust` = "#FC8D62", `MAD:hclust` = "#FC8D62", 
`CV:pam` = "#8DA0CB", `SD:pam` = "#8DA0CB", `MAD:pam` = "#8DA0CB", 
`CV:mclust` = "#E78AC3", `SD:mclust` = "#E78AC3", `MAD:mclust` = "#E78AC3", 
`CV:NMF` = "#A6D854", `SD:NMF` = "#A6D854", `MAD:NMF` = "#A6D854", 
`CV:skmeans` = "#FFD92F", `SD:skmeans` = "#FFD92F", `MAD:skmeans` = "#FFD92F", 
`CV:kmeans` = "#E5C494", `SD:kmeans` = "#E5C494", `MAD:kmeans` = "#E5C494"
)

library(dendextend)
library(ComplexHeatmap)
figure_d = grid.grabExpr(gridGraphics::grid.echo(function() {
	par(mar = c(4, 2, 2, 7), cex = 0.9)
	if(project == "recount2") {
		dend = rotate(dend, names(label_color))
	}
	color_labels(dend, col = label_color[labels(dend)]) %>% plot(horiz = TRUE, main = "Similarity between partitions")

	if(project == "GDS") {
		rect(0-0.05, 1-0.3, dend_heights(dend[[1]]) + 0.2, 6 + 0.3, border = label_color[1], lty = 2)
		rect(0-0.05, 7-0.3, dend_heights(dend[[2]][[1]]) + 0.2, 9 + 0.3, border = label_color[7], lty = 2)
		rect(0-0.05, 10-0.3, dend_heights(dend[[2]][[2]][[1]]) + 0.2, 12 + 0.3, border = label_color[10], lty = 2)
		rect(0-0.05, 13-0.3, dend_heights(dend[[2]][[2]][[2]][[1]]) + 0.2, 15 + 0.3, border = label_color[13], lty = 2)
		rect(0-0.05, 16-0.3, dend_heights(dend[[2]][[2]][[2]][[2]][[1]]) + 0.2, 18 + 0.3, border = label_color[16], lty = 2)
		rect(0-0.05, 19-0.3, dend_heights(dend[[2]][[2]][[2]][[2]][[2]][[1]]) + 0.2, 21 + 0.3, border = label_color[19], lty = 2)
		rect(0-0.05, 22-0.3, dend_heights(dend[[2]][[2]][[2]][[2]][[2]][[2]]) + 0.2, 24 + 0.3, border = label_color[22], lty = 2)
	}
}))
	

dev.off()

if(0) {

project = "GDS"
e1 = new.env()
sys.source("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_cohort/compare_methods.R", envir = e1)
figure_a1 = e1$figure_a
figure_b1 = e1$figure_b
figure_c1 = e1$figure_c
figure_d1 = e1$figure_d

project = "recount2"
e2 = new.env()
sys.source("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_cohort/compare_methods.R", envir = e2)
figure_a2 = e2$figure_a
figure_b2 = e2$figure_b
figure_c2 = e2$figure_c
figure_d2 = e2$figure_d

# dl = untangle(e1$dend, e2$dend, method = "step2side")
# figure_d1 = ggdendrogram(dl[[1]], rotate = TRUE) + scale_y_continuous(expand = c(0.01, 0.01), labels = NULL) +
#     ggtitle("Similarity of Partitions")

# figure_d2 = ggdendrogram(dl[[2]], rotate = TRUE) + scale_y_continuous(expand = c(0.01, 0.01), labels = NULL) +
#     ggtitle("Similarity of Partitions")

library(cowplot)
theme_set(theme_grey())

pdf(qq("@{BASE_DIR}/image/figure4.pdf"), width = 20, height = 20/2)
p = plot_grid(
		plot_grid(figure_a1, figure_b1, figure_c1, figure_d1, labels = c("A", "B", "C", "D"), align = "h", nrow = 1, rel_widths = c(1, 1.2, 1, 1)),
    	plot_grid(figure_a2, figure_b2, figure_c2, figure_d2, labels = c("E", "F", "G", "H"), align = "h", nrow = 1, rel_widths = c(1, 1.2, 1, 1)),
    	align = "v", ncol = 1
    )
print(p)
dev.off()


### running time for each dataset

collect_running_time = function(project) {
	if(project == "recount2") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
	} else if(project == "GDS") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
	}

	df = file.info(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id_list}", collapse = FALSE))
	id_list = basename(rownames(df)[df$isdir])

	
	method = NULL
	running_time = NULL
	n_sample = NULL
	pid = NULL
	i = 0
	for(id in id_list) {
		i = i + 1
		qqcat("loading @{id} (@{i}/@{length(id_list)})...\n")
		oe = try({
			res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_cola_all.rds"))

			pm = gsub("^.*:", "", names(res@list))
			x = sapply(res@list, function(x) x@running_time)
			x = tapply(x, pm, mean)
			running_time = c(running_time, x)
			method = c(method, names(x))
			n_sample = c(n_sample, rep(ncol(res), length(x)))
			pid = c(pid, rep(id, length(x)))
		})
	}

	data.frame(partition_method = method,
		       n_sample = n_sample,
		       running_time = running_time,
		       pid = pid,
		       stringsAsFactors = FALSE)
}

df1 = collect_running_time("GDS"); df1$dataset = "GDS"
df2 = collect_running_time("recount2"); df2$dataset = "recount2"
df = rbind(df1, df2)

df$partition_method = factor(df$partition_method, levels = c("hclust", "kmeans", "pam", "mclust", "skmeans", "NMF"))

library(ggplot2)
library(scales)
pdf(qq("@{BASE_DIR}/image/datasets_running_time.pdf"), width = 10, height = 5)
p = ggplot(df, aes(x = n_sample, y = running_time, col = partition_method)) +
	geom_point() + geom_smooth() +
	scale_y_continuous(trans = log10_trans(),
    	labels = trans_format("log10", function(x) paste0(10^x, "s"))) +
	scale_x_continuous(breaks = seq(0, 500, 100)) +
	facet_grid(.~ dataset, scales = "free_x", space = "free_x") +
	labs(x = "Number of samples", y = "Running time", col = "Partition method")
print(p)
dev.off()

}

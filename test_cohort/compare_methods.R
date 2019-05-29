

library(cola)
library(GetoptLong)
BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"

### test the running time for each method

register_NMF()

n = 20
tm_df = NULL
for(n in c(20, 50, 100, 200)) {
	for(pm in all_partition_methods()) {
		for(k in 2:6) {
			tm = NULL
			cat(paste0(n, "-", pm, "-", k), "...\n")
			f = cola:::get_partition_method(pm)
			counter = set_counter(100)
			for(i in 1:100) {
				m = matrix(runif(n*n), nr = n)
				t1 = nanotime(Sys.time()) 
				f(m, k)
				t2 = nanotime(Sys.time()) 
				tm[i] = as.numeric(t2 - t1)
				counter()
			}

			tm_df = rbind(tm_df, data.frame(n = n, pm = pm, k = k, time = tm))
		}
	}
}
save(tm_df, file = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/partition_method_running_time.RData")

pdf(qq("@{BASE_DIR}/image/partition_method_running_time.pdf"), width = 14, height = 5)

tm_df$pm = factor(tm_df$pm, levels = c("hclust", "kmeans", "pam", "mclust", "skmeans", "NMF"))
library(ggplot2)
library(scales)
ggplot(tm_df, aes(x=pm, y=time, fill=paste0("k=",k))) + 
    geom_boxplot() +
    scale_y_continuous(trans = log10_trans(),
    	 labels = trans_format("log10", function(x) paste0(10^(x-9), "s"))) +
    facet_wrap(~n, nrow = 1)
dev.off()



collect_stat = function(project) {
	if(project == "recount2") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
	} else if(project == "GDS") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
	}
	
	lt = list()
	for(id in id_list) {
		qqcat("loading @{id}...\n")
		oe = try({
			res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_cola_all.rds"))
			tb = suggest_best_k(res_list)
			tb = tb[order(rownames(tb)), ]
		})

		if(!inherits(oe, "try-error")) {
			lt[[id]] = tb
		}
	}
	return(lt)
}


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(matrixStats)
make_plot = function(lt, field) {
	# anno = data.frame(
	# 	top_value_method = gsub("^(.*):.*$", "\\1", rownames(lt[[1]])),
	# 	partition_method = gsub("^.*:(.*)$", "\\1", rownames(lt[[1]])),
	# 	stringsAsFactors = FALSE
	# )
	# anno_col = list(
	# 	top_value_method = structure(brewer.pal(4, "Set1"), names = unique(anno$top_value_method)),
	# 	partition_method = structure(brewer.pal(6, "Set2"), names = unique(anno$partition_method))
	# )

	m = do.call(rbind, lapply(lt, function(x) x[, field]))
	if(field == "PAC") {
		m = 1 - m
		field = "1 - PAC"
	}
	colnames(m) = rownames(lt[[1]])

	best_k = do.call(rbind, lapply(lt, function(x) x[, "best_k"]))
	colnames(best_k) = rownames(lt[[1]])
	p_best_k = t(apply(best_k, 2, function(x) {
		y = rep(0, 5)
		names(y) = 2:6
		tb = table(x)
		y[names(tb)] = tb
		y/length(x)
	}))

	ht = densityHeatmap(m, column_order = order(colMedians(m, na.rm = TRUE)), ylab = field, ylim = c(min(m, na.rm = TRUE), 1),
		# bottom_annotation = HeatmapAnnotation(df = anno, col = anno_col, gp = gpar(col = "white")),
		bottom_annotation = HeatmapAnnotation(p = anno_barplot(p_best_k, gp = gpar(fill = 2:6)), height = unit(3, "cm")),
		column_names_rot = 45,
		column_names_gp = gpar(col = ifelse(grepl("ATC|skmeans", colnames(m)), "red", "black")))
	draw(ht)

	return(invisible(m))
}

lt = collect_stat("GDS")


### compare best_k
m = matrix(NA, nrow = nrow(lt[[1]]), ncol = length(lt))
rownames(m) = rownames(lt[[1]])
colnames(m) = names(lt)
for(i in seq_along(lt)) {
	l = lt[[i]][, 6] != ""
	m[l, i] = lt[[i]][l, 1]
}
m2 = apply(m, 2, function(x) {
	x - median(x, na.rm = TRUE)
})
row_split = ifelse(grepl("ATC|skmeans", rownames(m)), "ATC|skmeans", "others")
# Heatmap(m, row_split = row_split, column_order = order(colSums(m)),
# 	cluster_rows = FALSE, cluster_columns = FALSE,
# 	col = colorRamp2(c(2, 6), c("white", "red")))
# dev.off()

x = colMeans(m[row_split == "ATC|skmeans", ], na.rm = T)
y = colMeans(m[row_split != "ATC|skmeans", ], na.rm = T)

p = t.test(x - y, alternative = "greater")$p.value

od = order(x - y)
pdf(qq("@{BASE_DIR}/image/recount2_best_k_diff.pdf"))
plot(sort(x -y), type= "p", main = p)
abline(h = 0)
abline(v = mean(which(sort(x - y) == 0)))
dev.off()

# cophcor
# PAC
# mean_silhouette
# concordance 
image1 = grid.grabExpr(make_plot(lt, "cophcor"))
image2 = grid.grabExpr(make_plot(lt, "PAC"))
image3 = grid.grabExpr(make_plot(lt, "mean_silhouette"))
image4 = grid.grabExpr(make_plot(lt, "concordance"))
pdf(qq("@{BASE_DIR}/image/GDS_stat_best_k_density_distribution.pdf"), width = 30, height = 8)
grid.newpage()
pushViewport(viewport(x = 0, width = 0.25, just = "left"))
grid.draw(image1)
popViewport()
pushViewport(viewport(x = 0.25, width = 0.25, just = "left"))
grid.draw(image2)
popViewport()
pushViewport(viewport(x = 0.5, width = 0.25, just = "left"))
grid.draw(image3)
popViewport()
pushViewport(viewport(x = 0.75, width = 0.25, just = "left"))
grid.draw(image4)
popViewport()
dev.off()


# compare different statistics
collect_all_stat = function(project) {
	if(project == "recount2") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
	} else if(project == "GDS") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
	}
	
	lt = list()
	for(id in id_list) {
		qqcat("loading @{id}...\n")
		oe = try({
			res = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_cola_all.rds"))

			df = do.call(rbind, lapply(2:6, function(k) {
				df = get_stats(res, k = k)
				df = as.data.frame(df)
				df$method = rownames(df)
				df$id = id
				best_k = guess_best_k(res)[rownames(df), ]
				df$best_k = best_k$best_k
				df$best_k_sig = best_k[, ncol(best_k)]
				rownames(df) = NULL
				df$PAC = 1 - df$PAC
				names(df)[names(df) == "PAC"] = "1-PAC"
				df$n_sample = ncol(res)
				df
			}))
		})

		if(!inherits(oe, "try-error")) {
			lt = c(lt, list(df))
		}
	}
	return(do.call("rbind", lt))
}

df_all = collect_all_stat("recount2")

library(GGally)
for(method in unique(df_all$method)) {
	l = df_all$method == method
	p = ggpairs(df_all[l, 2:5], mapping = aes(color = paste0("k=", df_all$k[l]))) + ggtitle(method)
	print(p)
}

foo = tapply(seq_len(nrow(df_all)), df_all[, c("k", "method")], function(ind) {
	cm = cor(df_all[ind, c("cophcor", "1-PAC", "mean_silhouette", "concordance")], method = "spearman")
	c("PAC_vs_cophcor" = cm["1-PAC", "cophcor"],
	  "mean_silhouette_vs_cophcor" = cm["mean_silhouette", "cophcor"],
	  "concordance_vs_cophcor" = cm["concordance", "cophcor"],
	  "mean_silhouette_vs_PAC" = cm["mean_silhouette", "1-PAC"],
	  "concordance_vs_PAC" = cm["concordance", "1-PAC"],
	  "concordance_vs_mean_silhouette" = cm["concordance", "mean_silhouette"])
})

arr = array(0, dim = c(nrow(foo), ncol(foo), length(foo[1, 1][[1]])),
	dimnames = list(rownames(foo), colnames(foo), names(foo[1, 1][[1]])))
for(i in seq_len(nrow(foo))) {
	for(j in seq_len(ncol(foo))) {
		arr[i, j, ] = foo[i, j][[1]]
	}
}

df = melt(arr)
colnames(df) = c("k", "method", "comparison", "correlation")
pdf(qq("@{BASE_DIR}/image/GDS_stat_correlation.pdf"), width = 20, height = 6.5)

ggplot(df, aes(x = comparison, y = correlation)) + 
	geom_boxplot() + facet_wrap( ~ factor(k), nrow = 1) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



## simply distribution of statistics
df = rbind(data.frame(k = df_all$k, value = df_all$cophcor, method = df_all$method, stat = "cophcor", n_sample = df_all$n_sample),
	       data.frame(k = df_all$k, value = df_all[, "1-PAC"], method = df_all$method, stat = "1-PAC", n_sample = df_all$n_sample),
	       data.frame(k = df_all$k, value = df_all$mean_silhouette, method = df_all$method, stat = "mean_silhouette", n_sample = df_all$n_sample),
	       data.frame(k = df_all$k, value = df_all$concordance, method = df_all$method, stat = "concordance", n_sample = df_all$n_sample)
)
pdf(qq("@{BASE_DIR}/image/recount2_stat_distribution.pdf"), width = 16, height = 9)

ggplot(df, aes(x = value, group = stat, fill = stat)) +
	geom_density() +
	facet_grid(stat ~ k)
ggplot(df, aes(y = value, x = factor(round(n_sample*2, -2)/2))) +
	geom_boxplot() +
	facet_grid(stat ~ k)
dev.off()


#### go through consensus heatmap and member ship heatmap to determine a proper PAC cutoff
make_cc_plot = function(project) {
	if(project == "recount2") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
	} else if(project == "GDS") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
	}
	par(ask = TRUE)
	on.exit(par(ask = FALSE))

	lt = list()
	for(id in id_list) {
		qqcat("loading @{id}...\n")
		oe = try({
			res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_cola_all.rds"))
		})

		if(!inherits(oe, "try-error")) {
			for(i in seq_along(res_list@list)) {
				stat = get_stats(res_list@list[[i]], k = guess_best_k(res_list@list[[i]]))
				consensus_heatmap(res_list@list[[i]], k = guess_best_k(res_list@list[[i]]))
				decorate_heatmap_body("Consensus", {
					grid.text(qq("pac = @{stat[1, 'PAC']}\nmean_silhouette = @{stat[1, 'mean_silhouette']}"),
						x = 0, y = 0, just = c("left", "bottom"))
				})
			}
		}
	}
	return(lt)
}
make_cc_plot("recount2")
make_cc_plot("GDS")


### check the functional ...

collect_p = function(project, field = c("gene", "BP", "CC", "MF")) {
	if(project == "recount2") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
	} else if(project == "GDS") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
	}

	field = match.arg(field)[1]
	
	lt = list()
	for(id in id_list) {
		qqcat("loading @{id}...\n")
		oe = try({
			load(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_sig_gene_GO.RData"))
			p = get(qq("p_sig_@{field}"))
			x = p[, 1]/p[, 2]
		})

		if(!inherits(oe, "try-error")) {
			lt[[id]] = x
		}
	}
	do.call(rbind, lt)
}
image = list()
for(field in c("gene", "BP", "MF", "CC")) {
	m = collect_p("GDS", field)
	ht = densityHeatmap(m, column_order = order(colMeans(m, na.rm = TRUE)), ylab = NULL, ylim = range(m, na.rm = TRUE),
		column_names_gp = gpar(col = ifelse(grepl("ATC|skmeans", colnames(m)), "red", "black")), show_heatmap_legend = FALSE,
		column_title = qq("density heatmap of sig_@{field}%"))
	image[[field]] = grid.grabExpr(draw(ht))
}
pdf(qq("@{BASE_DIR}/image/GDS_p_gene_GO.pdf"), width = 18.5, height = 5.2)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 4)))
for(i in seq_along(image)) {
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
	grid.draw(image[[i]])
	popViewport()
}
popViewport()
dev.off()

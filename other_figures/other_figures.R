library(cola)
library(GetoptLong)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(RColorBrewer)

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"

############ figure 2 ###############
rl = readRDS("/desktop-home/guz/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup.rds")

all_top_value_list = rl@.env$all_top_value_list[rl@top_value_method]

m = get_matrix(rl)

col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))

library(cowplot)
pdf(qq("@{BASE_DIR}/image/HSMM_single_cell_figure2.pdf"), width = 15, height = 10)

pl = list()
i = 0
for(nm in names(all_top_value_list)) {
	i = i + 1
	qqcat("making heatmap for @{nm}\n")
	od = order(all_top_value_list[[nm]], decreasing = TRUE)[1:2000]
	m2 = m[od, ]
	m2 = t(scale(t(m2)))
	
	ht_opt$TITLE_PADDING = unit(2, "mm")
	ht = Heatmap(m2, col = col_fun, column_title = qq("Top 2000 rows by @{nm}"),
		show_row_dend = FALSE, show_column_dend = FALSE, show_row_names = FALSE, 
		show_column_names = FALSE, show_heatmap_legend = FALSE, use_raster = TRUE)
	pl[[i]] = grid.grabExpr(draw(ht, newpage = FALSE, padding = unit(c(2, 2, 2, 2), "mm")))
}

lt = lapply(all_top_value_list, function(x) order(x, decreasing = TRUE)[1:2000])
col = RColorBrewer::brewer.pal(4, "Set2")
pl[[5]] = grid.grabExpr(grid.draw(plot(eulerr::euler(lt), fills = col)))

lgd = packLegend(
	Legend(col_fun = col_fun, title = "scaled expression"),
	Legend(title = "top-value method",
		labels = names(all_top_value_list), legend_gp = gpar(fill = col))
)
pl[[6]] = grid.grabExpr(draw(lgd, x = unit(0, "npc") + unit(5, "mm"), just = "left"))


plot_grid(pl[[1]], pl[[2]], pl[[3]], pl[[4]], pl[[5]], pl[[6]],
	nrow = 2, labels = LETTERS[1:6])
dev.off()



############ figure 3 ##############
rl = readRDS("/desktop-home/guz/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup.rds")

plot_ecdf = function(object, ...) {
	par(mar = c(4.1, 4.1, 4.1, 6), xpd = NA)
	plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "consensus value (x)", ylab = "P(X <= x)")
	for(i in seq_along(object@k)) {
		consensus_mat = get_consensus(object, k = object@k[i])
		f = ecdf(consensus_mat[lower.tri(consensus_mat)])
		x = seq(0, 1, length = 100)
		y = f(x)
		x = c(0, x)
		y = c(0, y)
		lines(x, y, col = i)
	}
	title("ECDF of consensus matrix")
	legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = paste0("k = ", object@k), 
					lty = 1, lwd = 2,
					col = seq_along(object@k), xjust = 0, yjust = 0.5,
					bty = "n")
}

plot_stat = function(x) {
	par(mar = c(4.1, 4.1, 4.1, 6), xpd = NA)
	plot(2:6, x, type = "b", xlab = "k (number of classes)", ylab = "1-PAC");
	points(3, x[2], col = "red", pch = 16)
	title("Determine the best number of classes")
	legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = "The best k", 
					pch = 16,
					col = "red", xjust = 0, yjust = 0.5,
					bty = "n")
}

library(cowplot)
p1 = grid.grabExpr(collect_plots(rl, k = 3, fun = consensus_heatmap))
ov = cola:::STAT_USED
assignInNamespace("STAT_USED", "1-PAC", ns = "cola")
p2 = grid.grabExpr(collect_stats(rl, k = 3, layout_nrow = 1))
assignInNamespace("STAT_USED", ov, ns = "cola")
p3 = grid.grabExpr(top_rows_overlap(rl, top_n = 500, method = "euler"))
res = rl["ATC:skmeans"]
p4 = ~plot_ecdf(res)
p5 = ~{plot_stat(get_stats(res)[, "1-PAC"]); abline(h = 0.9, lty = 2, col = "grey")}
p6 = grid.grabExpr(membership_heatmap(res, k = 3,))
p7 = ~dimension_reduction(res, k = 3, method = "PCA")
p8 = grid.grabExpr({set.seed(123); get_signatures(res, k = 3, simplify = TRUE)})
p9 = grid.grabExpr(collect_classes(rl, k = 3, show_row_dend = FALSE, merge_legends = TRUE))

# p = plot_grid(p2, p3, labels = c("B", "C"), ncol = 1)
# p = plot_grid(p1, p, labels = c("A", ""), nrow = 1, rel_widths = c(2, 1))
# p = plot_grid(p,
# 	plot_grid(p4, p5, p6, labels = c("D", "E", "F"), nrow = 1),
# 	plot_grid(p7, p8, p9, labels = c("G", "H", "I"), nrow = 1),
# 	ncol = 1, rel_heights = c(2 ,1.2, 1.2))

p = plot_grid(p1, plot_grid(p2, p9, p3, p4, p5, p6, p7, p8, nrow = 2, labels = LETTERS[2:9]),
	labels = c("A", ""), nrow = 2, rel_heights = c(2.5, 2))

pdf(qq("@{BASE_DIR}/image/Golub_leukemia_figure6.pdf"), width = 24, height = 24)
print(p)
dev.off()


########## figure 7 ########

rh = readRDS("/desktop-home/guz/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_hierarchical_partition.rds")
if(file.exists("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/HSMM_single_cell_ATC_skmeans.rds")) {
	res = readRDS("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/HSMM_single_cell_ATC_skmeans.rds")
} else {
	set.seed(123)
	res = consensus_partition(get_matrix(rh), max_k = 8, top_value_method = "ATC", partition_method = "skmeans", mc.cores = 4)
	saveRDS(res, file = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/HSMM_single_cell_ATC_skmeans.rds")
}

collect_classes(rh, anno = data.frame(CC = factor(get_classes(res, k = suggest_best_k(res))[, 1])))
dimension_reduction(rh)
get_signatures(rh)

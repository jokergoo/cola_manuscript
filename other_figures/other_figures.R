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

anno = get_anno(rl)
anno = anno[, c("Hours", "Media")]
anno_col = get_anno_col(rl)[c("Hours", "Media")]

pl = list()
i = 0
for(nm in names(all_top_value_list)) {
	i = i + 1
	qqcat("making heatmap for @{nm}\n")
	od = order(all_top_value_list[[nm]], decreasing = TRUE)[1:1000]
	m2 = m[od, ]
	m2 = t(scale(t(m2)))
	
	ht_opt$TITLE_PADDING = unit(2, "mm")
	ht = Heatmap(m2, col = col_fun, column_title = qq("Top 1000 rows by @{nm}"),
		top_annotation = HeatmapAnnotation(df = anno, col = anno_col, show_legend = FALSE, 
			show_annotation_name = FALSE, simple_anno_size = unit(3, "mm")),
		show_row_dend = FALSE, show_column_dend = FALSE, show_row_names = FALSE, 
		show_column_names = FALSE, show_heatmap_legend = FALSE, use_raster = TRUE)
	pl[[i]] = grid.grabExpr(draw(ht, newpage = FALSE, padding = unit(c(2, 2, 2, 2), "mm")))
}

lt = lapply(all_top_value_list, function(x) order(x, decreasing = TRUE)[1:1000])
col = RColorBrewer::brewer.pal(4, "Set2")
pl[[5]] = grid.grabExpr(grid.draw({
	gb = plot(eulerr::euler(lt), fills = col)
	gb$vp$width = unit(1, "npc") - unit(8, "mm")
	gb$vp$height = unit(1, "npc") - unit(8, "mm")
	print(gb)
}))

lgd = packLegend(
	Legend(title = "Hours", labels = names(anno_col$Hours), direction = "horizontal", nrow = 1, legend_gp = gpar(fill = anno_col$Hours)),
	Legend(title = "Media", labels = names(anno_col$Media), direction = "horizontal", nrow = 1, legend_gp = gpar(fill = anno_col$Media)),
	Legend(col_fun = col_fun, title = "Scaled expression", direction = "horizontal"),
	Legend(title = "Top-value method", labels = names(all_top_value_list), direction = "horizontal", nrow = 2, legend_gp = gpar(fill = col)),
	direction = "horizontal", max_width = unit(14/4, "inch")
)
pl[[6]] = grid.grabExpr(draw(lgd, x = unit(0, "npc") + unit(5, "mm"), just = "left"))


id_mapping = map_to_entrez_id("ENSEMBL")
lt = rl@.env$all_top_value_list

gl = lapply(lt, function(x) {
	gene = rownames(rl)[order(x, decreasing = TRUE)[1:1000]]
	fl = functional_enrichment(gene, id_mapping = function(x) id_mapping[gsub("\\.\\d+$", "", x)])$BP
	fl$ID[fl$p.adjust < 0.01]
})

library(simplifyEnrichment)

go_id = unique(unlist(gl))
mm = matrix(0, ncol = 4, nrow = length(go_id), dimnames = list(go_id, names(gl)))
mm[go_id %in% gl$SD, 1] = 1
mm[go_id %in% gl$CV, 2] = 1
mm[go_id %in% gl$MAD, 3] = 1
mm[go_id %in% gl$ATC, 4] = 1
pl[[6]] = grid.grabExpr(simplifyGO(GO_similarity(go_id, "BP"), 
	word_cloud_grob_param = list(max_width = unit(130, "mm")),
	exclude_words = c("process", "regulation"), control = list(cutoff = 0.97),
	fontsize_range = c(5, 16),
	ht_list = Heatmap(mm, col = c("0" = NA, "1" = "darkgreen"), width = unit(16, "mm"), 
		cluster_columns = FALSE, show_heatmap_legend = FALSE)),
    width = 10.5, height = 8/2.2*1.2)

pdf(qq("@{BASE_DIR}/image/HSMM_single_cell_figure2.pdf"), width = 14, height = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 4, heights = c(1, 1.2))))
for(i in 1:4) {
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
	grid.draw(pl[[i]])
	popViewport()
}
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
pushViewport(viewport(y = 0.9, just = "top", height = 0.6))
grid.text("Overlap of top 1000 genes", y = unit(1, "npc") - unit(1.25, "mm"), just = "top",
	gp = gpar(fontsize = 14))
grid.draw(pl[[5]])
popViewport()
draw(lgd, y = unit(0.05, "npc"), just = "bottom")
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
pushViewport(viewport(x = 0, width = unit(10.5, "inch"), height = unit(8/2.2*1.2, "inch"), just = "left"))
grid.draw(pl[[6]])
popViewport(3)
dev.off()


### compare euler plot and upset plot
pdf(qq("@{BASE_DIR}/image/HSMM_single_cell_figure2_euler_upset.pdf"), width = 10, height = 5)

pl = list()
lt = lapply(all_top_value_list, function(x) order(x, decreasing = TRUE)[1:1000])
col = RColorBrewer::brewer.pal(4, "Set2")
pl[[1]] = grid.grabExpr(grid.draw({
	gb = plot(eulerr::euler(lt), fills = col)
	gb$vp$width = unit(1, "npc") - unit(8, "mm")
	gb$vp$height = unit(1, "npc") - unit(8, "mm")
	print(gb)
}))

cm = make_comb_mat(lt)
pl[[2]] = grid.grabExpr(draw(UpSet(cm)))


plot_grid(pl[[1]], pl[[2]], nrow = 1, rel_widths = c(1, 1.5))
dev.off()

############ figure 3 ##############
rl = readRDS("/icgc/dkfzlsdf/analysis/B080/guz/repo/cola_examples/Golub_leukemia/Golub_leukemia_subgroup.rds")

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
p1 = grid.grabExpr(collect_plots(rl, k = 3, fun = consensus_heatmap, simplify = TRUE))
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
p8 = grid.grabExpr({set.seed(88); get_signatures(res, k = 3, simplify = TRUE)})
p9 = grid.grabExpr(collect_classes(rl, k = 3, show_row_dend = FALSE, merge_legends = TRUE, simplify = TRUE))

# p = plot_grid(p2, p3, labels = c("B", "C"), ncol = 1)
# p = plot_grid(p1, p, labels = c("A", ""), nrow = 1, rel_widths = c(2, 1))
# p = plot_grid(p,
# 	plot_grid(p4, p5, p6, labels = c("D", "E", "F"), nrow = 1),
# 	plot_grid(p7, p8, p9, labels = c("G", "H", "I"), nrow = 1),
# 	ncol = 1, rel_heights = c(2 ,1.2, 1.2))


pdf(qq("@{BASE_DIR}/image/Golub_leukemia_figure6.pdf"), width = 18, height = 18*1.323529)
p = plot_grid(p1, plot_grid(p2, p9, p3, p4, p5, p6, p7, p8, nrow = 3, labels = LETTERS[2:9]),
	labels = c("A", ""), nrow = 2, rel_heights = c(1.5, 2))
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

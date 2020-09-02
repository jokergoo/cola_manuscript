library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)
library(GetoptLong)
library(ggplot2)
library(cola)

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19", 
	package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
probe = IlluminaHumanMethylation450kanno.ilmn12.hg19 # change to a short name

library(GEOquery)
if(file.exists("/desktop-home/guz/project/development/cola_examples/GBM_450K/GSE36278_450K.RData")) {
	load("/desktop-home/guz/project/development/cola_examples/GBM_450K/GSE36278_450K.RData")
} else {
	gset = getGEO("GSE36278")
	save(gset, file = "/desktop-home/guz/project/development/cola_examples/GBM_450K/GSE36278_450K.RData")
}

mat = exprs(gset[[1]])
colnames(mat) = phenoData(gset[[1]])@data$title
mat = mat[rownames(getAnnotation(probe, what = "Locations")), ]

l = getAnnotation(probe, what = "Locations")$chr %in% paste0("chr", 1:22) & 
	is.na(getAnnotation(probe, what = "SNPs.137CommonSingle")$Probe_rs)
mat = mat[l, ]

cpg_anno = getAnnotation(probe, "Islands.UCSC")$Relation_to_Island[l]
cpg_anno2 = cpg_anno
cpg_anno2[ cpg_anno == "Island" ] = "island"
cpg_anno2[ cpg_anno %in% c("N_Shelf", "S_Shelf", "OpenSea") ] = "sea"
cpg_anno2[ cpg_anno %in% c("N_Shore", "S_Shore") ] = "shore"

cgi_col = c("island" = "red", "shore" = "blue", "sea" = "grey")
tb = table(cpg_anno2)[c("island", "shore", "sea")]

mat1 = as.matrix(mat[, grep("GBM", colnames(mat))])   # tumor samples
colnames(mat1) = gsub("GBM", "dkfz", colnames(mat1))

mat1[is.na(mat1)] = runif(sum(is.na(mat1)))

phenotype = read.table("/desktop-home/guz/project/development/cola_examples/GBM_450K/450K_annotation.txt", header = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
phenotype = phenotype[colnames(mat1), ]
colnames(phenotype)[1:2] = c("dkfz_subtype", "tcga_subtype")

anno_col = list(
    dkfz_subtype = structure(names = c("IDH", "K27", "G34", "RTK I PDGFRA", "Mesenchymal", "RTK II Classic"), brewer.pal(6, "Set1"))
)

if(!file.exists("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GBM_450k_top_value_list_sd_ATC.rds")) {
	library(bsub)
	bsub_chunk({
		top_value_list = list()
		top_value_list$SD = cola:::get_top_value_method("SD")(mat1)
		top_value_list$ATC = cola:::get_top_value_method("ATC")(mat1, mc.cores = 8)
	    saveRDS(top_value_list, file = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GBM_450k_top_value_list_sd_ATC.rds")
	}, name = "meth_top_value_list", variables = "mat1", hour = 10, memory = 10, core = 8, packages = "cola")
}

top_value_list = readRDS("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GBM_450k_top_value_list_sd_ATC.rds")

meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
i = 1
od = order(top_value_list[[i]], decreasing = TRUE)[1:5000]
ht = Heatmap(mat1[od, ], name = "meth", column_title = qq("top @{length(od)} rows of @{names(top_value_list)[i]}"),
    show_row_names = FALSE, show_row_dend = FALSE, show_column_names = FALSE,
    clustering_method_rows = "ward.D2",
    row_split = cpg_anno2[od], col = meth_col_fun, cluster_row_slices = FALSE,
    left_annotation = rowAnnotation(cpg_anno = cpg_anno2[od], col = list(cpg_anno = cgi_col)),
    bottom_annotation = HeatmapAnnotation(dkfz_subtype = phenotype$dkfz_subtype, col = anno_col))
figure_a = grid.grabExpr(draw(ht, merge_legend = TRUE))
i = 2
od = order(top_value_list[[i]], decreasing = TRUE)[1:5000]
ht = Heatmap(mat1[od, ], name = "meth", column_title = qq("top @{length(od)} rows of @{names(top_value_list)[i]}"),
    show_row_names = FALSE, show_row_dend = FALSE, show_column_names = FALSE,
    clustering_method_rows = "ward.D2",
    row_split = cpg_anno2[od], col = meth_col_fun, cluster_row_slices = FALSE,
    left_annotation = rowAnnotation(cpg_anno = cpg_anno2[od], col = list(cpg_anno = cgi_col)),
    bottom_annotation = HeatmapAnnotation(dkfz_subtype = phenotype$dkfz_subtype, col = anno_col))
figure_b = grid.grabExpr(draw(ht, merge_legend = TRUE))


######### figure c ###############

rl_all = readRDS("/desktop-home/guz/project/development/cola_examples/GBM_450K/GBM_450K_cgi_all_subgroup.rds")
rl_island = readRDS("/desktop-home/guz/project/development/cola_examples/GBM_450K/GBM_450K_cgi_island_subgroup.rds")
rl_shore = readRDS("/desktop-home/guz/project/development/cola_examples/GBM_450K/GBM_450K_cgi_shore_subgroup.rds")
rl_sea = readRDS("/desktop-home/guz/project/development/cola_examples/GBM_450K/GBM_450K_cgi_sea_subgroup.rds")

method = "SD:skmeans"
class_df = data.frame(
    all = get_classes(rl_all[method], k = 6)[, 1],
    all_sil = get_classes(rl_all[method], k = 6)[, 3],
    island = get_classes(rl_island[method], k = 6)[, 1],
    island_sil = get_classes(rl_island[method], k = 6)[, 3],
    shore = get_classes(rl_shore[method], k = 6)[, 1],
    shore_sil = get_classes(rl_shore[method], k = 6)[, 3],
    sea = get_classes(rl_sea[method], k = 4)[, 1],
    sea_sil = get_classes(rl_sea[method], k = 4)[, 3]
)

ref_class = class_df$all
class_df$all = relabel_class(class_df$all, ref_class, return_map = FALSE)
map_island = relabel_class(class_df$island, ref_class)
class_df$island = as.numeric(map_island[as.character(class_df$island)])
map_shore = relabel_class(class_df$shore, ref_class)
class_df$shore = as.numeric(map_shore[as.character(class_df$shore)])
map_sea = relabel_class(class_df$sea, ref_class); map_sea = map_sea[1:4]
class_df$sea = as.numeric(map_sea[as.character(class_df$sea)])

adjust_by_transparency = function(col, transparency) {
    rgb( 1 - (1 - t(col2rgb(col)/255)) * (1 - transparency))
}

switch_labels = function(x, label1, label2) {
	x[x == label2] = 0
	x[x == label1] = label2
	x[x == 0] = label1
	x
}

class_df$all = switch_labels(class_df$all, 2, 5)
class_df$island = switch_labels(class_df$island, 2, 5)
class_df$shore = switch_labels(class_df$shore, 2, 5)
class_df$sea = switch_labels(class_df$sea, 2, 5)

class_df$all = switch_labels(class_df$all, 5, 6)
class_df$island = switch_labels(class_df$island, 5, 6)
class_df$shore = switch_labels(class_df$shore, 5, 6)
class_df$sea = switch_labels(class_df$sea, 5, 6)

subtype_order = c("G34" = 1, "IDH" = 2, "K27" = 3,
	              "Mesenchymal" = 6, "RTK I PDGFRA" = 5, "RTK II Classic" = 4)

class_mat = t(as.matrix(class_df[, c(1, 3, 5, 7)]))
silhouette_mat = t(as.matrix(class_df[, c(1, 3, 5, 7) + 1]))
silhouette_mat[silhouette_mat < 0] = 0

ht = Heatmap(class_mat[-1, ], name = "Class",
    col = cola:::brewer_pal_set2_col[1:7], cluster_rows = FALSE,
    column_order = order(class_df$island, class_df$sea, class_df$shore, subtype_order[phenotype$dkfz_subtype]),
    bottom_annotation = HeatmapAnnotation(dkfz_subtype = phenotype$dkfz_subtype, 
        col = anno_col, show_legend = TRUE),
    column_title = qq("partitions from @{method}"),
    rect_gp = gpar(type = "none"),
    row_title = NULL,
    layer_fun = function(j, i, x, y, w, h, fill) {
        col = fill
        # col = adjust_by_transparency(fill, 1 - pindex(silhouette_mat, j, i))
        grid.rect(x, y, w, h, gp = gpar(fill = col, col = col))
    })
figure_c = grid.grabExpr(draw(ht))

### figure D

method = "ATC:skmeans"
class_df = data.frame(
    all = get_classes(rl_all[method], k = 3)[, 1],
    all_sil = get_classes(rl_all[method], k = 3)[, 3],
    island = get_classes(rl_island[method], k = 5)[, 1],
    island_sil = get_classes(rl_all[method], k = 5)[, 3],
    shore = get_classes(rl_shore[method], k = 3)[, 1],
    shore_sil = get_classes(rl_all[method], k = 3)[, 3],
    sea = get_classes(rl_sea[method], k = 4)[, 1],
    sea_sil = get_classes(rl_all[method], k = 4)[, 3]
)

ref_class = class_df$all
class_df$all = relabel_class(class_df$all, ref_class, return_map = FALSE)
map_island = relabel_class(class_df$island, ref_class)
class_df$island = as.numeric(map_island[as.character(class_df$island)])
map_shore = relabel_class(class_df$shore, ref_class)
class_df$shore = as.numeric(map_shore[as.character(class_df$shore)])
map_sea = relabel_class(class_df$sea, ref_class); map_sea = map_sea[1:4]
class_df$sea = as.numeric(map_sea[as.character(class_df$sea)])

class_df$island = switch_labels(class_df$island, 4, 1)

class_mat = t(as.matrix(class_df[, c(1, 3, 5, 7)]))
silhouette_mat = t(as.matrix(class_df[, c(1, 3, 5, 7) + 1]))
silhouette_mat[silhouette_mat < 0] = 0

ht = Heatmap(class_mat[-1, ], name = "Class",
    col = cola:::brewer_pal_set2_col[1:7], cluster_rows = FALSE,
    column_order = order(class_df$all, class_df$sea, class_df$shore, subtype_order[phenotype$dkfz_subtype]),
    bottom_annotation = HeatmapAnnotation(dkfz_subtype = phenotype$dkfz_subtype, 
        col = anno_col, show_legend = TRUE),
    column_title = qq("partitions from @{method}"),
    rect_gp = gpar(type = "none"),
    row_title = NULL,
    layer_fun = function(j, i, x, y, w, h, fill) {
        col = fill
        # col = adjust_by_transparency(fill, 1 - pindex(silhouette_mat, j, i))
        grid.rect(x, y, w, h, gp = gpar(fill = col, col = col))
    })
figure_d = grid.grabExpr(draw(ht))


library(cowplot)
p = plot_grid(plot_grid(figure_a, figure_b, nrow = 1, labels = c("A", "B")), 
    figure_c, figure_d,
	labels = c("", "C", "D"), nrow = 3, rel_heights = c(1, 0.25, 0.25))

pdf(qq("@{BASE_DIR}/image/methylation_figure7.pdf"), width = 10, height = 9)
print(p)
dev.off()


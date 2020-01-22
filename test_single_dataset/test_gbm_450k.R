options(showWarnCalls = TRUE, showErrorCalls = TRUE)

library(GetoptLong)

cgi = "all"
GetoptLong("cgi=s", "island|shore|sea")


setwd("/desktop-home/guz/project/development/cola_examples/GBM_450K/")

library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)

library(cola)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19", 
	package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
probe = IlluminaHumanMethylation450kanno.ilmn12.hg19 # change to a short name

library(GEOquery)
if(file.exists("GSE36278_450K.RData")) {
	load("GSE36278_450K.RData")
} else {
	gset = getGEO("GSE36278")
	save(gset, file = "GSE36278_450K.RData")
}

mat = exprs(gset[[1]])
colnames(mat) = phenoData(gset[[1]])@data$title
mat = mat[rownames(getAnnotation(probe, what = "Locations")), ]

l = getAnnotation(probe, what = "Locations")$chr %in% paste0("chr", 1:22) & 
	is.na(getAnnotation(probe, what = "SNPs.137CommonSingle")$Probe_rs)
mat = mat[l, ]

cgi_anno = getAnnotation(probe, "Islands.UCSC")$Relation_to_Island[l]
table(cgi_anno)

mat1 = as.matrix(mat[, grep("GBM", colnames(mat))])   # tumor samples
colnames(mat1) = gsub("GBM", "dkfz", colnames(mat1))

phenotype = read.table("450K_annotation.txt", header = TRUE, sep = "\t", row.names = 1, 
	check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
phenotype = phenotype[colnames(mat1), 1:2]
colnames(phenotype) = c("dkfz_subtype", "tcga_subtype")

mat1[is.na(mat1)] = runif(sum(is.na(mat1)))

cgi_regexp = cgi
if(cgi == "all") {
	cgi_regexp = ".*"
} else if(cgi %in% c("shelf", "sea")) {
	cgi_regexp = "shelf|sea"
}
mat1 = mat1[grepl(cgi_regexp, cgi_anno, ignore.case = TRUE), , drop = FALSE]


anno_col = list(
	dkfz_subtype = structure(names = c("IDH", "K27", "G34", "RTK I PDGFRA", "Mesenchymal", "RTK II Classic"), brewer.pal(6, "Set1")),
    tcga_subtype = structure(names = c("G-CIMP+", "Cluster #2", "Cluster #3"), brewer.pal(3, "Set1"))
)

set.seed(123)
rl = run_all_consensus_partition_methods(
	mat1, 
	top_n = seq(min(2000, round(nrow(data) * 0.1)), min(10000, round(nrow(data) * 0.5)), length.out = 5),
	max_k = 10,
	scale_rows = FALSE, 
	anno = phenotype, 
	anno_col = anno_col, 
	mc.cores = 4)

saveRDS(rl, file = qq("GBM_450K_cgi_@{cgi}_subgroup.rds"))

cola_opt(group_diff = 0.1)
cola_report(rl, output_dir = qq("GBM_450K_cgi_@{cgi}_subgroup_cola_report"), mc.cores = 4)

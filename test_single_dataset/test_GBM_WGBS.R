
anno_text =
"id    subtype
AK015  IDH
AK041  IDH
AK066  IDH
AK068  IDH
AK076  IDH
AK085  IDH
AK102  IDH
AK103  IDH
AK124  IDH
AK199  IDH
AK213  IDH
AK231  IDH
AK005  MES
AK006  MES
AK030  MES
AK055  MES
AK071  MES
AK072  MES
AK079  MES
AK081  MES
AK088  MES
AK091  MES
AK134  MES
AK139  MES
AK153  MES
AK185  MES
AK188  MES
AK195  MES
AK218  MES
AK227  MES
AK236  MES
AK002  RTK_I
AK003  RTK_I
AK043  RTK_I
AK049  RTK_I
AK051  RTK_I
AK142  RTK_I
AK149  RTK_I
AK156  RTK_I
AK165  RTK_I
AK173  RTK_I
AK183  RTK_I
AK217  RTK_I
AK035  RTK_II
AK053  RTK_II
AK074  RTK_II
AK089  RTK_II
AK098  RTK_II
AK099  RTK_II
AK100  RTK_II
AK117  RTK_II
AK123  RTK_II
AK132  RTK_II
AK133  RTK_II
AK158  RTK_II
AK167  RTK_II
AK178  RTK_II
AK205  RTK_II
AK216  RTK_II
AK226  RTK_II
"

SAMPLE = read.table(textConnection(anno_text), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
SAMPLE_ID = rownames(SAMPLE)
SUBTYPE_COLOR = RColorBrewer::brewer.pal(4, "Set1")
names(SUBTYPE_COLOR) = c("IDH", "MES", "RTK_I", "RTK_II")
COLOR = list(subtype = SUBTYPE_COLOR)

library(GetoptLong)
library(GenomicRanges)
library(data.table)
library(epik)  # https://github.com/jokergoo/epik

methylation_hooks$get_by_chr = function(chr) {
	qqcat("reading /icgc/dkfzlsdf/analysis/B080/guz/cola_test/GSE121721_GBM/GSE121721_CG_methylation_bsseq_smoothed_@{chr}.txt.gz\n")
    tb = fread(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GSE121721_GBM/GSE121721_CG_methylation_bsseq_smoothed_@{chr}.txt.gz"), header = TRUE)
    tb = as.data.frame(tb)[, SAMPLE_ID]
    gr = fread(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GSE121721_GBM/GSE121721_CG_positions_@{chr}.txt.gz"), header = TRUE)
    gr = as.data.frame(gr)
    gr = GRanges(seqnames = gr[, 1], ranges = IRanges(gr[, 2], gr[, 3]))

    obj2 = list(gr = gr, meth = tb)
    return(obj2)
}

CHROMOSOME = paste0("chr", 1:22)
GENOME = "hg19"

CGI = read.table(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GSE121721_GBM/cpgIslandExt.bed"), stringsAsFactors = FALSE)
CGI = CGI[CGI[, 1] %in% CHROMOSOME, ]
CGI = GRanges(seqnames = CGI[[1]], ranges = IRanges(CGI[[2]], CGI[[3]]))
CGI_SHORE = setdiff(flank(CGI, 2000, both = TRUE), CGI)

GENOMIC_FEATURE_LIST = list()
GENOMIC_FEATURE_LIST$cgi = CGI
GENOMIC_FEATURE_LIST$shore = CGI_SHORE

library(GenomicFeatures)
chrom_info = getChromInfoFromUCSC(GENOME)
chrom_info = chrom_info[chrom_info[, 1] %in% CHROMOSOME, ]
gr_chrom = GRanges(seqnames = chrom_info[, 1], ranges = IRanges(1, chrom_info[, 2]))

GENOMIC_FEATURE_LIST$sea = setdiff(gr_chrom, c(CGI, CGI_SHORE))
GENOMIC_FEATURE_LIST$all = gr_chrom

windowsize = 1000
library(EnrichedHeatmap)
gr_list = lapply(GENOMIC_FEATURE_LIST, function(gr) {
	gr = makeWindows(gr, w = windowsize, short.keep = TRUE)
	gr[width(gr) >= windowsize/2]
})

set.seed(123)
gr_list = lapply(gr_list, function(gr) {
	if(length(gr) > 500000) {
		p = 500000/length(gr)
		gr = gr[sample(c(TRUE, FALSE), length(gr), p = c(p, 1-p), replace = TRUE)]
	}
	gr
})

gr_list = get_mean_methylation_in_genomic_features(SAMPLE_ID, chromosome = CHROMOSOME, genomic_features = gr_list)

mat_list = lapply(gr_list, function(gr) {
	m = as.matrix(mcols(gr))
	m = m[apply(m, 1, function(x) all(!is.na(x))), grep("mean_meth", colnames(m))]
	m
})

save(mat_list, SAMPLE, COLOR, file = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GSE121721_GBM/GSE121721_GBM_processed.RData")

### apply cola analysis ###
library(cola)
library(bsub) # https://github.com/jokergoo/bsub
for(nm in names(mat_list)) {

	bsub_chunk(name = qq("cola_GSE121721_GBM_@{nm}"), hour = 50, memory = 30, core = 4, variable = c("nm"), 
		packages = c("cola", "GetoptLong"), 
		image = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GSE121721_GBM/GSE121721_GBM_processed.RData",
	{
		set.seed(123)
		mat = mat_list[[nm]]
	
		rl = run_all_consensus_partition_methods(
			mat, 
			top_n = seq(min(2000, round(nrow(mat) * 0.1)), min(10000, round(nrow(mat) * 0.5)), length.out = 5),
			max_k = 10,
			scale_rows = FALSE,
			anno = SAMPLE, 
			anno_col = COLOR,
			mc.cores = 4
		)
		saveRDS(rl, file = qq("/desktop-home/guz/project/development/cola_examples/GSE121721_GBM_WGBS/GSE121721_GBM_WGBS_@{nm}_subgroup.rds"))
		cola_report(rl, output_dir = qq("/desktop-home/guz/project/development/cola_examples/GSE121721_GBM_WGBS/GSE121721_GBM_WGBS_@{nm}_subgroup_cola_report"))
	})
}

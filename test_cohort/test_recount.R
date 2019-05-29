library(GetoptLong)

rerun = FALSE
GetoptLong(
	"pid=s", "pid",
	"rerun!", "rerun"
)


######################################################
## get the expression matrix
######################################################

library(SummarizedExperiment)

load(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_rse_gene.Rdata"))

count = assays(rse_gene)$counts
df = rowData(rse_gene)

if(ncol(count) > 400) {
	set.seed(123)
	ind = sample(ncol(count), 400)
	count = count[, ind]
}

######################################################
## TPM normalization
######################################################

tpm_normalize = function(count, gene_length) {
	all_count = colSums(count)

	tpm = matrix(0, nrow = nrow(count), ncol = ncol(count))
	rownames(tpm) = rownames(count)
	colnames(tpm) = colnames(count)
	for(i in seq_len(nrow(count))) {
		tpm[i, ] = count[i, ] / gene_length[i]
	}

	for(j in seq_len(ncol(count))) {
		tpm[, j] = tpm[, j] / all_count[j]
	}

	tpm = 10^9 * tpm
	return(tpm)
}
tpm = tpm_normalize(count, df$bp_length)

load("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gencode_v25_transcript_merged.RData")

### only the protein coding genes
gt = sapply(gene_annotation$gtf, function(x) x$type)
l = gt[rownames(count)] == "protein_coding"

count = count[l, , drop = FALSE]
tpm = tpm[l, , drop = FALSE]

## remove genes that are not expressed
l = apply(count, 1, function(x) sum(x > 2)/length(x) > 0.5)

mat = log2(tpm[l, , drop = FALSE] + 1)

##############################################
## apply cola analysis
##############################################

library(cola)
register_NMF()
mat = adjust_matrix(mat)

if(!rerun && file.exists(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_all.rds"))) {
 	res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_all.rds"))
} else {
	if(file.exists(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_all.rds"))) {
		file.remove(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_all.rds"))
	}

	if(nrow(mat) > 6000) {
		res_list = run_all_consensus_partition_methods(mat, mc.cores = 4, top_n = c(1000, 2000, 3000))
	} else {
		res_list = run_all_consensus_partition_methods(mat, mc.cores = 4)
	}
	saveRDS(res_list, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_all.rds"))
}
cola_report(res_list, output_dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_report"), mc.cores = 4)

#################################################
## check the proportion of significant genes and GO terms
#################################################

load("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/id_mapping_entrez_id_hash.RData")
library(clusterProfiler) # entrez id as central id

cutoff = 0.05
p_sig_gene = matrix(NA, nr = length(res_list@list), nc = 2)
p_sig_BP = matrix(NA, nr = length(res_list@list), nc = 2)
p_sig_MF = matrix(NA, nr = length(res_list@list), nc = 2)
p_sig_CC = matrix(NA, nr = length(res_list@list), nc = 2)
rownames(p_sig_gene) = rownames(p_sig_BP) = rownames(p_sig_MF) = rownames(p_sig_CC) = names(res_list@list)
colnames(p_sig_gene) = colnames(p_sig_BP) = colnames(p_sig_MF) = colnames(p_sig_CC) = c("n_sig", "n_all")

for(i in seq_along(res_list@list)) {
	cat("=====================================\n")
	qqcat("* signature genes and GO terms for @{names(res_list@list)[i]}\n")
	sig_df = get_signatures(res_list@list[[i]], k = suggest_best_k(res_list@list[[i]]), plot = FALSE, verbose = FALSE)
	if(is.null(sig_df)) {
		sig_gene = NULL
	} else {
		p_sig_gene[i, ] = c(sum(sig_df$fdr < cutoff), nrow(get_matrix(res_list)))
		sig_gene = rownames(sig_df[sig_df$fdr < cutoff, , drop = FALSE])
	}

	if(length(sig_gene)) {
		sig_gene = gsub("\\.\\d+$", "", sig_gene)
		sig_gene = ENSEMBL2ENTREZID[sig_gene]
		sig_gene = sig_gene[!is.na(sig_gene)]
		sig_gene = unique(sig_gene)

		if(length(sig_gene)) {
			cat("  - gene set enrichment, GO:BP\n")
			ego = enrichGO(gene = sig_gene,
                OrgDb         = 'org.Hs.eg.db',
                ont           = "BP",
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 1000,
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
        		readable      = TRUE)
			ego = as.data.frame(ego)
			p_sig_BP[i, ] = c(sum(ego$p.adjust < cutoff), nrow(ego))

			cat("  - gene set enrichment, GO:MF\n")
			ego = enrichGO(gene = sig_gene,
                OrgDb         = 'org.Hs.eg.db',
                ont           = "MF",
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 1000,
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
        		readable      = TRUE)
			ego = as.data.frame(ego)
			p_sig_MF[i, ] = c(sum(ego$p.adjust < cutoff), nrow(ego))

			cat("  - gene set enrichment, GO:CC\n")
			ego = enrichGO(gene = sig_gene,
                OrgDb         = 'org.Hs.eg.db',
                ont           = "CC",
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 1000,
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
        		readable      = TRUE)
			ego = as.data.frame(ego)
			p_sig_CC[i, ] = c(sum(ego$p.adjust < cutoff), nrow(ego))
		}
	}
}

save(p_sig_gene, p_sig_BP, p_sig_MF, p_sig_CC, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_sig_gene_GO.RData"))

# for(pid in dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/")) {
# 	cmd = qq("module load R/3.3.1; Rscript /icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_cohort/test_recount.R --pid @{pid}")
# 	dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/")
# 	cmd = qq("perl /desktop-home/guz/project/development/ngspipeline2/bsub_single_line.pl --hour 50 --memory 20 --core 4 --name cola_recount_@{pid} --dir @{dir} --command '@{cmd}' --enforce")
# 	system(cmd)
# }

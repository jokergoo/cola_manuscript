library(GetoptLong)

cores = 4
GetoptLong(
	"pid=s", "pid",
	"cores=i", "cores"
)


######################################################
## get the expression matrix
######################################################

library(SummarizedExperiment)

load(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_rse_gene.Rdata"))

count = assays(rse_gene)$counts
df = rowData(rse_gene)

if(ncol(count) > 500) {
	set.seed(123)
	ind = sample(ncol(count), 500)
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

res_list = run_all_consensus_partition_methods(mat, mc.cores = cores)
saveRDS(res_list, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_all.rds"))
cola_report(res_list, output_dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_report"), 
	title = qq("cola Report for recount2:@{pid}"), mc.cores = cores)


library(GetoptLong)

rerun = FALSE
GetoptLong(
	"id=s", "gds id",
	"rerun!", "rerun"
)

library(GEOquery)


######################################################
## get the expression matrix
######################################################
gds = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}.rds"))

ind = which(is.na(gds@dataTable@table[, 1]))
if(length(ind)) {
	qqcat("find @{length(ind)} rows have NA IDs.\n")
	gds@dataTable@table = gds@dataTable@table[-ind, ]
}
oe = try(eset <- GDS2eSet(gds))
if(inherits(oe, "try-error")) eset = GDS2eSet(gds, getGPL=FALSE)
mat = exprs(eset)

######################################################
## test whether need to perform log2 transformation
######################################################

## if in `value_type` column contain "log", most probably it means it is a two-channel
## microarray and the values are log fold change
if(!grepl("log", gds@header$value_type)) {
	mat[mat < 0] = 0
	x = mat[mat >= 0]
	ks_stat = ks.test(x, "pnorm", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))$statistic

	if(ks_stat > 0.3) {
		mat = log2(mat + 1)
	}
}

if(id == "GDS2808") {
	mat = mat[, apply(mat, 2, function(x) sum(is.na(x))/length(x) < 0.25), drop = FALSE]
}

library(cola)
register_NMF()
mat = adjust_matrix(mat)


################################################
## apply quantile normalization 
################################################
if(!grepl("log", gds@header$value_type)) {
	library(preprocessCore)
	cn = colnames(mat)
	rn = rownames(mat)
	mat = normalize.quantiles(mat)
	colnames(mat) = cn
	rownames(mat) = rn
}


##############################################
## apply cola analysis
##############################################
if(!rerun && file.exists(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_all.rds"))) {
	res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_all.rds"))
} else {
	if(file.exists(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_all.rds"))) {
		file.remove(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_all.rds"))
	}

	if(nrow(mat) > 6000) {
		res_list = run_all_consensus_partition_methods(mat, mc.cores = 4, top_n = c(1000, 2000, 3000))
	} else {
		res_list = run_all_consensus_partition_methods(mat, mc.cores = 4)
	}
	saveRDS(res_list, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_all.rds"))
}
cola_report(res_list, output_dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_report"), mc.cores = 4)


#################################################
## check the proportion of significant genes and GO terms
#################################################

library(clusterProfiler) # entrez id as central id

cutoff = 0.05
p_sig_gene = matrix(NA, nr = length(res_list@list), nc = 2)
p_sig_BP = matrix(NA, nr = length(res_list@list), nc = 2)
p_sig_MF = matrix(NA, nr = length(res_list@list), nc = 2)
p_sig_CC = matrix(NA, nr = length(res_list@list), nc = 2)
rownames(p_sig_gene) = rownames(p_sig_BP) = rownames(p_sig_MF) = rownames(p_sig_CC) = names(res_list@list)
colnames(p_sig_gene) = colnames(p_sig_BP) = colnames(p_sig_MF) = colnames(p_sig_CC) = c("n_sig", "n_all")

platform = gds@header$platform
id_mapping = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/@{platform}_probe_id_to_entrez_id.rds"))


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
		sig_gene = id_mapping[sig_gene]
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

save(p_sig_gene, p_sig_BP, p_sig_MF, p_sig_CC, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}//@{id}_sig_gene_GO.RData"))


# for(id in dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/", pattern = "GDS")) {
# 	cmd = qq("module load R/3.3.1; Rscript /icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_cohort/test_GDS.R --id @{id}")
# 	dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}")
# 	cmd = qq("perl /desktop-home/guz/project/development/ngspipeline2/bsub_single_line.pl --hour 50 --memory 20 --core 4 --name cola_@{id} --dir @{dir} --command '@{cmd}' --enforce")
# 	system(cmd)
# }

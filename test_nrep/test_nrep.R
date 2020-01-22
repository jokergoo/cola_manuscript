library(GetoptLong)

dataset = "TCGA_GBM"
by = "row"
GetoptLong(
	"run=i", "run",
	"rep=i", "rep",
	"dataset=s", "dataset",
	"by=s", "row|column"
)

options(showWarnCalls = TRUE, showErrorCalls = TRUE)

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"
root = "/desktop-home/guz"

library(cola)

if(dataset == "TCGA_GBM") {
	m = read.table(qq("@{root}/project/development/cola_examples/TCGA_GBM/unifiedScaled.txt"), header = TRUE, row.names = 1, check.names = FALSE)
	m = as.matrix(m)

	subtype = read.table(qq("@{root}/project/development/cola_examples/TCGA_GBM/TCGA_unified_CORE_ClaNC840.txt"), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])

	m = m[, names(subtype)]
	m = adjust_matrix(m)

} else if(dataset == "HSMM_single_cell") {
	library(HSMMSingleCell)

	data(HSMM_expr_matrix)

	# `HSMM_expr_matrix` is a FPKM matrix
	m = adjust_matrix(log10(HSMM_expr_matrix + 1))
	
	gt = readRDS(qq("@{root}/project/development/cola_examples/HSMM_single_cell/gene_type_gencode_v17.rds"))
	m = m[gt[rownames(m)] == "protein_coding", , drop = FALSE]
}

register_NMF()

rl = run_all_consensus_partition_methods(
	m, 
	partition_repeat = rep,
	mc.cores = 4,
	sample_by = by
)
if(by == "column") dataset = paste0(dataset, "_by_column")
dir.create(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}/"), showWarnings = FALSE)
saveRDS(rl, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}/@{dataset}_subgroup_repeat_@{rep}_run_@{run}.rds"))

# for(dataset in c("TCGA_GBM", "HSMM_single_cell")) {
# 	for(rep in c(25, 50, 100, 200)) {
# 		for(run in 1:100) {
# 			by = "column"
# 			if(file.exists(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}_by_column/@{dataset}_by_column_subgroup_repeat_@{rep}_run_@{run}.rds"))) next
# 			cmd = qq("module load R/3.6.0; Rscript /icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_nrep/test_nrep.R --rep @{rep} --run @{run} --dataset @{dataset} --by @{by}")
# 			dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/test_partition_repeats/@{dataset}_by_column")
# 			cmd = qq("perl /desktop-home/guz/project/development/ngspipeline2/bsub_single_line.pl --hour 100 --memory 5 --core 4 --name @{dataset}_subgroup_repeat_@{rep}_run_@{run} --dir @{dir} --command '@{cmd}' --enforce")
# 			system(cmd)
# 		}
# 	}
# }



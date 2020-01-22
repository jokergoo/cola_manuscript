
library(GetoptLong)
library(cola)

dataset = "TCGA_GBM"
k = 1
GetoptLong(
	"dataset=s", "dataset",
	"k=i", k
)

rl = readRDS(qq("/desktop-home/guz/project/development/cola_examples/@{dataset}/@{dataset}_subgroup.rds"))
m = get_matrix(rl)

dir.create(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}/"), showWarnings = FALSE)

rl = run_all_consensus_partition_methods(m, top_value_method = c("SD", "ATC"),
	partition_method = c("hclust", "skmeans"), sample_by = "row")
saveRDS(rl, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}/@{dataset}_resampling_by_row_@{k}.rds"))

rl = run_all_consensus_partition_methods(m, top_value_method = c("SD", "ATC"),
	partition_method = c("hclust", "skmeans"), sample_by = "column")
saveRDS(rl, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}/@{dataset}_resampling_by_column_@{k}.rds"))


# library(bsub)
# for(dataset in c("Golub_leukemia", "HSMM_single_cell", "MCF10CA_scRNAseq", "Ritz_ALL", "TCGA_GBM")) {
# 	bsub_opt$output_dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/row_or_column_resampling/@{dataset}")
# 	for(k in 1:100) {
# 		bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_row_column_resampling/row_or_column_resampling.R",
# 			name = qq("@{dataset}_row_or_column_resampling_@{k}"), argv = qq("--dataset @{dataset} --k @{k}"),
# 			hour = 10, memory = 10)
# 	}
# }



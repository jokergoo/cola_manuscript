options(showWarnCalls = TRUE, showErrorCalls = TRUE)

setwd("/desktop-home/guz/project/development/cola_examples/HSMM_single_cell/")

library(cola)
library(RColorBrewer)

library(HSMMSingleCell)

data(HSMM_expr_matrix)
data(HSMM_sample_sheet)

# `HSMM_expr_matrix` is a FPKM matrix
m = adjust_matrix(log2(HSMM_expr_matrix + 1))
anno = HSMM_sample_sheet[, c("Hours", "Media", "State")]
anno_col = list(
	Hours = structure(brewer.pal(9, "Blues")[c(2, 4, 6, 8)], names = c("0", "24", "48", "72")),
	Media = c("GM" = "orange", "DM" = "purple"),
	State = c("1" = "red", "2" = "blue", "3" = "green"))

gt = readRDS("gene_type_gencode_v17.rds")
m = m[gt[rownames(m)] == "protein_coding", , drop = FALSE]

register_NMF()

set.seed(123)
rl = run_all_consensus_partition_methods(
	m,
	mc.cores = 4, 
	anno = anno,
	anno_col = anno_col
)

saveRDS(rl, file = "HSMM_single_cell_subgroup.rds")
cola_report(rl, output_dir = "HSMM_single_cell_subgroup_cola_report", mc.cores = 4)

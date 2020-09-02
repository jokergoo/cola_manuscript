options(showWarnCalls = TRUE, showErrorCalls = TRUE)

setwd("/desktop-home/guz/project/development/cola_examples/MCF10CA_scRNAseq/")

library(cola)

tpm = readRDS("MCF10CA_scRNAseq_tpm.rds")

m = log2(tpm + 1)

cell_type = ifelse(grepl("round", colnames(m)), "round", "aberrant")
cell_col = c("aberrant" = "red", "round" = "blue")

m = adjust_matrix(m)

register_NMF()

set.seed(123)
rl = run_all_consensus_partition_methods(
	m, 
	mc.cores = 4,
	anno = data.frame(cell_type = cell_type), 
	anno_col = list(cell_type = cell_col)
)

saveRDS(rl, file = "MCF10CA_scRNAseq_subgroup.rds")
cola_report(rl, output_dir = "MCF10CA_scRNAseq_subgroup_cola_report", mc.cores = 4)


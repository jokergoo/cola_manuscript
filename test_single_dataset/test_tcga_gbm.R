options(showWarnCalls = TRUE, showErrorCalls = TRUE)

setwd("/desktop-home/guz/project/development/cola_examples/TCGA_GBM/")

library(cola)
library(RColorBrewer)

m = read.table("unifiedScaled.txt", header = TRUE, row.names = 1, check.names = FALSE)
m = as.matrix(m)

subtype = read.table("TCGA_unified_CORE_ClaNC840.txt", sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])
subtype_col = structure(seq_len(4), names = unique(subtype))

m = m[, names(subtype)]
m = adjust_matrix(m)

register_NMF()

set.seed(123)
rl = run_all_consensus_partition_methods(
	m, 
	mc.cores = 4,
	anno = data.frame(subtype = subtype), 
	anno_col = list(subtype = subtype_col)
)
saveRDS(rl, file = "TCGA_GBM_subgroup.rds")
cola_report(rl, output_dir = "TCGA_GBM_subgroup_cola_report", mc.cores = 4)

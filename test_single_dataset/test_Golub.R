options(showWarnCalls = TRUE, showErrorCalls = TRUE)

setwd("/desktop-home/guz/project/development/cola_examples/Golub_leukemia/")

library(cola)

library(golubEsets)
data(Golub_Merge)
m = exprs(Golub_Merge)
colnames(m) = paste0("sample_", colnames(m))
anno = pData(Golub_Merge)

m[m <= 1] = NA
m = log10(m)

m = adjust_matrix(m)

library(preprocessCore)
cn = colnames(m)
rn = rownames(m)
m = normalize.quantiles(m)
colnames(m) = cn
rownames(m) = rn

register_NMF()

anno = anno[, "ALL.AML", drop = FALSE]
anno_col = list("ALL.AML" = c("ALL" = "red", "AML" = "blue"))

set.seed(123)
rl = run_all_consensus_partition_methods(
	m,
	mc.cores = 4, 
	anno = anno,
	anno_col = anno_col
)

saveRDS(rl, file = "Golub_leukemia_subgroup.rds")
cola_report(rl, output_dir = "Golub_leukemia_subgroup_cola_report", mc.cores = 4)

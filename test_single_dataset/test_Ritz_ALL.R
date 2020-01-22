options(showWarnCalls = TRUE, showErrorCalls = TRUE)

setwd("/desktop-home/guz/project/development/cola_examples/Ritz_ALL/")

library(cola)

library(ALL)
data(ALL)

m = exprs(ALL)
anno = pData(ALL)
anno = anno[, c("sex", "age", "BT")]

m = adjust_matrix(m)

library(preprocessCore)
cn = colnames(m)
rn = rownames(m)
m = normalize.quantiles(m)
colnames(m) = cn
rownames(m) = rn


register_NMF()


set.seed(123)
rl = run_all_consensus_partition_methods(
	m,
	mc.cores = 4, 
	anno = anno
)
saveRDS(rl, file = "Ritz_ALL_subgroup.rds")
cola_report(rl, output_dir = "Ritz_ALL_subgroup_cola_report", mc.cores = 4)

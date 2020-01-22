

library(cola)
library(GetoptLong)
BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"

project = "GDS"
GetoptLong("project=s", "project")

if(project == "recount2") {
	id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
} else if(project == "GDS") {
	id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
}

df = file.info(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id_list}", collapse = FALSE))
id_list = basename(rownames(df)[df$isdir])

commonness = function(lt) {
	n = length(lt)
	m = matrix(NA, nrow = n, ncol = n)
	for(i in 1:(n-1)) {
		for(j in (i+1):n) {
			m[i, j] = m[j, i] = commonness_single(lt[[i]], lt[[j]])
		}
	}
	m
}

commonness_single = function(x1, x2) {
	s1 = intersect(x1, x2)
	s2 = union(x1, x2)
	length(s1)/length(s2)
}

library(parallel)
df_list = mclapply(id_list, function(id) {	
	df = data.frame(id = character(0), method = character(0), k = character(0), commonness = numeric(0))	
	rl = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_cola_all.rds"))

	for(j in seq_along(rl@list)) {
		best_k = suggest_best_k(rl@list[[j]])
		if(!is.null(attr(best_k, "optional"))) {
			qqcat("####### @{id} @{names(rl@list)[j]} #######\n")
			all_best_k = c(best_k, attr(best_k, "optional"))
			sig_lt = list()
			for(ki in seq_along(all_best_k)) {
				sig_lt[[ ki ]] = get_signatures(rl@list[[j]], k = all_best_k[ki], plot = FALSE)[, 1]
				cat("\n")
			}
			df = rbind(df, data.frame(id = id, method = names(rl@list)[j], 
				commonness = max(commonness(sig_lt), na.rm = TRUE), 
				k = paste(all_best_k, collapse = ","),
				stringsAsFactors = FALSE))
		}
	}
	df
}, mc.cores = 4)

df = do.call("rbind", df_list)

saveRDS(df, file = qq("@{BASE_DIR}/@{project}/@{project}_compare_best_k.rds"))


# library(bsub)
# bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_cohort/compare_best_k.R",
# 	.project = "GDS", name = "cola_compare_best_k_GDS", memory = 10, hour = 20, core = 4)
# bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_cohort/compare_best_k.R",
# 	.project = "recount2", name = "cola_compare_best_k_recount2", memory = 10, hour = 20, core = 4)


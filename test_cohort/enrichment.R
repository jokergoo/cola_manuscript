
library(GetoptLong)

library(bsub)

#################################################
## GDS: the proportion of significant genes and GO terms
#################################################

for(id in dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/", pattern = "GDS\\d+")) {
	if(!file.info(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}"))$isdir) next
	bsub_opt$output_dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}")
	bsub_chunk(name = qq("GDS_@{id}_sig_gene_GO"), hour = 40, mem = 10, variables = "id", {

		library(GEOquery)
		library(cola)
		library(GetoptLong)
		
		gds = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}.rds"))

		platform = gds@header$platform
		id_mapping_vec = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/@{platform}_probe_id_to_entrez_id.rds"))

		id_mapping_fun = function(sig_gene) {
			sig_gene = id_mapping_vec[sig_gene]
			sig_gene = sig_gene[!is.na(sig_gene)]
			sig_gene = unique(sig_gene)
			sig_gene = unique(unlist(strsplit(sig_gene, " /// ")))
			sig_gene
		}

		res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_all.rds"))

		p_sig_gene = list()
		lt_GO = list()

		for(k in 2:6) {
			p_sig_gene[[as.character(k)]] = list()
			lt_GO[[as.character(k)]] = list()

			for(nm in names(res_list@list)) {
				res = res_list@list[[nm]]
				sig_df = get_signatures(res, k = k, plot = FALSE, verbose = FALSE)
				p_sig_gene[[as.character(k)]][[nm]] = nrow(sig_df)/nrow(res)
				lt_GO[[as.character(k)]][[nm]] = functional_enrichment(res, k = k, id_mapping = id_mapping_fun)
			}
		}

		save(p_sig_gene, lt_GO, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_sig_gene_GO.RData"))
	})
}


#################################################
## recount2: the proportion of significant genes and GO terms
#################################################
for(pid in dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/")) {
	if(!file.info(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}"))$isdir) next
	bsub_opt$output_dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/recount2/@{pid}")
	bsub_chunk(name = qq("recount_@{pid}_sig_gene_GO"), hour = 40, mem = 10, variables = "pid", {

		library(cola)
		library(GetoptLong)

		p_sig_gene = list()
		lt_GO = list()

		res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_all.rds"))

		for(k in 2:6) {
			p_sig_gene[[as.character(k)]] = list()
			lt_GO[[as.character(k)]] = list()

			for(nm in names(res_list@list)) {
				res = res_list@list[[nm]]
				sig_df = get_signatures(res, k = k, plot = FALSE, verbose = FALSE)
				p_sig_gene[[as.character(k)]][[nm]] = nrow(sig_df)/nrow(res)
				lt_GO[[as.character(k)]][[nm]] = functional_enrichment(res, k = k)
			}
		}

		save(p_sig_gene, lt_GO, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_sig_gene_GO.RData"))
	})
}


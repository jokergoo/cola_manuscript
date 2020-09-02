
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



collect_enrichment = function(project) {
	if(project == "recount2") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
	} else if(project == "GDS") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
	}

	df = file.info(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id_list}", collapse = FALSE))
	id_list = basename(rownames(df)[df$isdir])

	get_p_sig_gene_tb = function(p_sig_gene) {
		data.frame(k = unlist(lapply(names(p_sig_gene), function(x) rep(x, length(p_sig_gene[[x]])))),
			       method = unlist(lapply(p_sig_gene, names)),
			       p_sig = unlist(lapply(p_sig_gene, function(x) unlist(lapply(x, function(y) {
			       	if(length(y) == 0) {
			       		NA
			       	} else {
			       		y
			       	}
			       })))), stringsAsFactors = FALSE)
	}

	get_p_sig_go_tb = function(lt_GO) {
		data.frame(k = unlist(lapply(names(lt_GO), function(x) rep(x, length(lt_GO[[x]])))),
			       method = unlist(lapply(lt_GO, names)),
			       p_sig = unlist(lapply(lt_GO, function(x) unlist(lapply(x, function(y) {
			       		if(is.data.frame(y)) {
			       			sum(y$p.adjust < 0.05)/nrow(y)
			       		} else {
				       	    mean(sapply(y, function(z) {
				       	    	sum(z$p.adjust < 0.05)/nrow(z)
				       	    }), na.rm = TRUE)
				       	}
			       })))), stringsAsFactors = FALSE)
	}

	p_sig_gene_tb = data.frame(k = numeric(0), method = character(0), p_sig = numeric(0), id = character(0))
	p_sig_go_tb = data.frame(k = numeric(0), method = character(0), p_sig = numeric(0), id = character(0))
	i = 0
	for(id in id_list) {
		i = i + 1
		qqcat("loading @{id} (@{i}/@{length(id_list)})...\n")
		oe = try({
			load(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_sig_gene_GO.RData"))
			tb1 = get_p_sig_gene_tb(p_sig_gene)
			tb2 = get_p_sig_go_tb(lt_GO)
			tb1$pid = id
			tb2$pid = id

			p_sig_gene_tb = rbind(p_sig_gene_tb, tb1)
			p_sig_go_tb = rbind(p_sig_go_tb, tb2)
		})
	}

	list(p_sig_gene_tb = p_sig_gene_tb, p_sig_go_tb = p_sig_go_tb)
}

lt = collect_enrichment("GDS")
saveRDS(lt, file = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/go_list_GDS.rds")

lt = collect_enrichment("recount2")
saveRDS(lt, file = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/go_list_recount2.rds")


########
lt = readRDS("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/go_list_GDS.rds")
pdf("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/image/GDS_sig_gene_GO.pdf")
for(k in 2:6) {
	tb = lt$p_sig_gene_tb[lt$p_sig_gene_tb == k, ]
	foo = split(tb$p_sig, tb$method)
	ht = densityHeatmap(foo, column_order = order(sapply(foo, mean, na.rm =T)), ylim = c(0, 1),
		column_title = qq("Density heatmap of sig_gene%, k = @{k}"),
		column_names_gp = gpar(col = ifelse(grepl("skmeans|ATC", names(foo)), "red", "black")),
		ylab = "Percent of significant genes")
	draw(ht)

	tb = lt$p_sig_go_tb[lt$p_sig_go_tb == k, ]
	foo = split(tb$p_sig, tb$method)
	ht = densityHeatmap(foo, column_order = order(sapply(foo, mean, na.rm =T)), ylim = c(0, 1),
		column_title = qq("Density heatmap of sig_GO_BP%, k = @{k}"),
		column_names_gp = gpar(col = ifelse(grepl("skmeans|ATC", names(foo)), "red", "black")),
		ylab = "Percent of significant GO BP terms")
	draw(ht)
}
dev.off()

lt = readRDS("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/go_list_recount2.rds")
pdf("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/image/recount2_sig_gene_GO.pdf")
for(k in 2:6) {
	tb = lt$p_sig_gene_tb[lt$p_sig_gene_tb == k, ]
	foo = split(tb$p_sig, tb$method)
	ht = densityHeatmap(foo, column_order = order(sapply(foo, mean, na.rm =T)), ylim = c(0, 1),
		column_title = qq("Density heatmap of sig_gene%, k = @{k}"),
		column_names_gp = gpar(col = ifelse(grepl("skmeans|ATC", names(foo)), "red", "black")),
		ylab = "Percent of significant genes")
	draw(ht)

	tb = lt$p_sig_go_tb[lt$p_sig_go_tb == k, ]
	foo = split(tb$p_sig, tb$method)
	ht = densityHeatmap(foo, column_order = order(sapply(foo, mean, na.rm =T)), ylim = c(0, 1),
		column_title = qq("Density heatmap of sig_GO_BP%, k = @{k}"),
		column_names_gp = gpar(col = ifelse(grepl("skmeans|ATC", names(foo)), "red", "black")),
		ylab = "Percent of significant GO BP terms")
	draw(ht)
}
dev.off()

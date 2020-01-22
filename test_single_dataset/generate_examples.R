
library(GetoptLong)
library(bsub)

bsub_opt$output_dir = "/desktop-home/guz/project/development/cola_examples/"

bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_single_dataset/test_Golub.R", 
    name = "test_Golub", memory = 10, hour = 100, core = 4)
bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_single_dataset/test_tcga_gbm.R", 
    name = "test_tcga_gbm", memory = 10, hour = 100, core = 4)
bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_single_dataset/test_Ritz_ALL.R", 
    name = "test_Ritz_ALL", memory = 10, hour = 100, core = 4)
bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_single_dataset/test_HSMM_single_cell.R",
    name = "test_HSMM_single_cell", memory = 10, hour = 100, core = 4)
bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_single_dataset/test_MCF10CA_scRNAseq.R",
    name = "test_MCF10CA_scRNAseq", memory = 10, hour = 100, core = 4)

bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_single_dataset/test_gbm_450k.R",
    name = "test_gbm_450k", memory = 40, hour = 50, core = 4)
bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_single_dataset/test_gbm_450k.R",
    name = "test_gbm_450k_cgi_island", .cgi = "island", memory = 40, hour = 50, core = 4)
bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_single_dataset/test_gbm_450k.R",
    name = "test_gbm_450k_cgi_shore", .cgi = "shore", memory = 40, hour = 50, core = 4)
bsub_script("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/test_single_dataset/test_gbm_450k.R",
    name = "test_gbm_450k_cgi_sea", .cgi = "sea", memory = 40, hour = 50, core = 4)


####### only generate reports
bsub_opt$output_dir = "/desktop-home/guz/project/development/cola_examples/"
bsub_opt$packages = c("cola", "GetoptLong")
for(dataset in c("Golub_leukemia", "HSMM_single_cell", "MCF10CA_scRNAseq", "Ritz_ALL", "TCGA_GBM")) {
    bsub_chunk({
        rl = readRDS(qq("/desktop-home/guz/project/development/cola_examples/@{dataset}/@{dataset}_subgroup.rds"))
        cola_report(rl, output_dir = qq("/desktop-home/guz/project/development/cola_examples/@{dataset}/@{dataset}_subgroup_cola_report"), 
            title = qq("cola Report for @{dataset}"), mc.cores = 4)
    }, name = qq("cola_report_@{dataset}"), hour = 10, memory = 20, core = 4, variable = c("dataset"))
}
for(cgi in c("all", "island", "shore", "sea")) {
    bsub_chunk({
        rl = readRDS(qq("/desktop-home/guz/project/development/cola_examples/GBM_450K/GBM_450K_cgi_@{cgi}_subgroup.rds"))
        cola_report(rl, output_dir = qq("/desktop-home/guz/project/development/cola_examples/GBM_450K/GBM_450K_cgi_@{cgi}_subgroup_cola_report"), 
            title = qq("cola Report for GBM_450K cgi_@{cgi}"), mc.cores = 4)
    }, name = qq("cola_report_GBM_450K_cgi_@{cgi}"), hour = 10, memory = 20, core = 4, variable = c("cgi"))
}


##################################################
###  function enrichment
bsub_opt$packages = "cola"
# Golub_leukemia: HU6800, hu6800.db
bsub_chunk(name = "cola_Golub_GO", hour = 10, memory = 4, {
    library(hu6800.db)
    x <- hu6800ENTREZID
    mapped_probes <- mappedkeys(x)
    id_mapping <- unlist(as.list(x[mapped_probes]))
    rl = readRDS("/desktop-home/guz/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup.rds")
    lt = functional_enrichment(rl, id_mapping = id_mapping)
    saveRDS(lt, file = "/desktop-home/guz/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_lt_GO.rds")
})

# HSMM single cell: Ensembl
bsub_chunk(name = "cola_HSMM_single_cell_GO", hour = 10, memory = 4, {
    id_mapping = map_to_entrez_id("ENSEMBL")
    rl = readRDS("/desktop-home/guz/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup.rds")
    lt = functional_enrichment(rl, id_mapping = function(x) id_mapping[gsub("\\.\\d+$", "", x)])
    saveRDS(lt, file = "/desktop-home/guz/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_lt_GO.rds")
})

# MCF10CA, Ensembl, no version numbers
bsub_chunk(name = "cola_MCF10CA_GO", hour = 10, memory = 4, {
    id_mapping = map_to_entrez_id("ENSEMBL")
    rl = readRDS("/desktop-home/guz/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup.rds")
    lt = functional_enrichment(rl, id_mapping = function(x) id_mapping[gsub("\\.\\d+$", "", x)])
    saveRDS(lt, file = "/desktop-home/guz/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup_lt_GO.rds")
})

# Ritz ALL: HGU95aV2, hgu95av2.db
bsub_chunk(name = "cola_Ritz_ALL_GO", hour = 10, memory = 4, {
    library(hgu95av2.db)
    x <- hgu95av2ENTREZID
    mapped_probes <- mappedkeys(x)
    id_mapping <- unlist(as.list(x[mapped_probes]))
    rl = readRDS("/desktop-home/guz/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup.rds")
    lt = functional_enrichment(rl, id_mapping = id_mapping)
    saveRDS(lt, file = "/desktop-home/guz/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_lt_GO.rds")
})

# TCGA GBM, gene symbol
bsub_chunk(name = "cola_TCGA_GBM_GO", hour = 10, memory = 4, {
    library(org.Hs.eg.db)
    id_mapping = map_to_entrez_id("SYMBOL")
    rl = readRDS("/desktop-home/guz/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup.rds")
    lt = functional_enrichment(rl, id_mapping = id_mapping)
    saveRDS(lt, file = "/desktop-home/guz/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup_lt_GO.rds")
})



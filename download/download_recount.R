library(GetoptLong)

recount_abstract = read.csv("recount_selection_2018-11-13_14_42_40.csv", header = TRUE)

# sample size in 50 ~ 500
anno = recount_abstract[recount_abstract[, 2] >= 50 & recount_abstract[, 2] < 500, ]

library(SummarizedExperiment)

for(pid in anno$accession) {
	dir.create(pid)
	url = paste0("http://duffel.rail.bio/recount/v2/", pid, "/rse_gene.Rdata")
	download.file(url, destfile = qq("@{pid}/@{pid}_rse_gene.Rdata"))
}


# GTEx datasets
df = read.table(textConnection(
"http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_adrenal_gland.Rdata adrenal_gland
# http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_bladder.Rdata bladder
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_adipose_tissue.Rdata adipose_tissue
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_blood_vessel.Rdata blood_vessel
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_blood.Rdata blood
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_bone_marrow.Rdata bone_marrow
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_brain.Rdata brain
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_breast.Rdata breast
# http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_cervix_uteri.Rdata cervix_uteri
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_colon.Rdata colon
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_esophagus.Rdata esophagus
# http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_fallopian_tube.Rdata fallopian_tube
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_heart.Rdata heart
# http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_kidney.Rdata kidney
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_liver.Rdata liver
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_lung.Rdata lung
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_muscle.Rdata muscle
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_nerve.Rdata nerve
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_ovary.Rdata ovary
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_pancreas.Rdata pancreas
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_pituitary.Rdata pituitary
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_prostate.Rdata prostate
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_salivary_gland.Rdata salivary_gland
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_skin.Rdata skin
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_small_intestine.Rdata small_intestine
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_spleen.Rdata spleen
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_stomach.Rdata stomach
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_testis.Rdata testis
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_thyroid.Rdata thyroid
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_uterus.Rdata uterus
http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_vagina.Rdata vagina
"), stringsAsFactors = FALSE)

df[[2]] = paste0("GTEx_", gsub(" ", "_", df[[2]]))
for(i in seq_len(nrow(df))) {
	dir.create(df[i, 2])
	download.file(df[i, 1], destfile = qq("@{df[i, 2]}/@{df[i, 2]}_rse_gene.Rdata"))
}

# TCGA datasets
df = read.table(textConnection(
"http://duffel.rail.bio/recount/v2/TCGA/rse_gene_adrenal_gland.Rdata adrenal_gland
# http://duffel.rail.bio/recount/v2/TCGA/rse_gene_bile_duct.Rdata bile_duct
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_bladder.Rdata bladder
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_bone_marrow.Rdata bone_marrow
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_brain.Rdata brain
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_breast.Rdata breast
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_cervix.Rdata cervix
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_colorectal.Rdata colorectal
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_esophagus.Rdata esophagus
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_eye.Rdata eye
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_head_and_neck.Rdata head_and_neck
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_kidney.Rdata kidney
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_liver.Rdata liver
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_lung.Rdata lung
# http://duffel.rail.bio/recount/v2/TCGA/rse_gene_lymph_nodes.Rdata lymph_nodes
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_ovary.Rdata ovary
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_pancreas.Rdata pancreas
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_pleura.Rdata pleura
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_prostate.Rdata prostate
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_skin.Rdata skin
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_soft_tissue.Rdata soft_tissue
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_stomach.Rdata stomach
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_testis.Rdata testis
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_thymus.Rdata thymus
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_thyroid.Rdata thyroid
http://duffel.rail.bio/recount/v2/TCGA/rse_gene_uterus.Rdata uterus
"), stringsAsFactors = FALSE)

df[[2]] = paste0("TCGA_", gsub(" ", "_", df[[2]]))
for(i in seq_len(nrow(df))) {
	dir.create(df[i, 2])
	download.file(df[i, 1], destfile = qq("@{df[i, 2]}/@{df[i, 2]}_rse_gene.Rdata"))
}

library(GEOquery)
library(GetoptLong)

## get GPL information of GDS datasets
gpl = NULL
for(id in dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/", pattern = "GDS")) {
	gds = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}.rds"))
	gpl[id] = gds@header$platform
}

# download GPL data
for(id in unique(gpl)) {
	try({
		gds = getGEO(id, destdir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/"))
		saveRDS(gds, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/@{id}.rds"))
	})
}

##############################################################
### save mapping from ACCNUM/UNIGENE/SYMBOL/ENSEMBL/REFSEQ to entrez IDs
##############################################################
library(hash)
global_mapping_to_entrez_id = function(from) {
	x <- getFromNamespace(qq("org.Hs.eg@{from}"), ns = "org.Hs.eg.db")
	mapped_genes <- mappedkeys(x)
	xx <- as.list(x[mapped_genes])

	ENTREZID = rep(names(xx), times = sapply(xx, length))
	df = data.frame(ENTREZID, unlist(xx))
	names(df) = c("ENTREZID", from)
	df = df[sample(nrow(df), nrow(df)), ]
	df = df[!duplicated(df[, 2]), ]
	x = structure(df[, 1], names = df[, 2])
	return(x)
}

library(org.Hs.eg.db)
ACCNUM2ENTREZID = global_mapping_to_entrez_id("ACCNUM")
UNIGENE2ENTREZID = global_mapping_to_entrez_id("UNIGENE")
SYMBOL2ENTREZID = global_mapping_to_entrez_id("SYMBOL")
ENSEMBL2ENTREZID = global_mapping_to_entrez_id("ENSEMBL")
REFSEQ2ENTREZID = global_mapping_to_entrez_id("REFSEQ")

save(ACCNUM2ENTREZID, UNIGENE2ENTREZID, SYMBOL2ENTREZID, ENSEMBL2ENTREZID, REFSEQ2ENTREZID, 
	file = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/id_mapping_entrez_id_hash.RData")


##############################################################
## save mapping from probe IDs to entrez IDs
## If there is no direct mapping, probe IDs are first mapped to ACCNUM/UNIGENE/SYMBOL/ENSEMBL/REFSEQ, then to entrez id
## If there is multiple mapping, take a random one.
##############################################################

library(hashmap)
library(epik)
## a named vector from probe id to entrez id
probe_id_2_entrez_id = function(platform) {
	gpl = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/@{platform}.rds"))

	all_nm = gpl@dataTable@columns[, 1]
	if(any(grepl("(Entrez.*Gene.*ID|Entrez.*?Gene)", all_nm, ignore.case = TRUE))) {
		qqcat("@{platform}: ID -> ENTREZID\n")
		ind = grep("(Entrez.*Gene.*ID|Entrez.*?Gene)", all_nm, ignore.case = TRUE)
		if(any(grepl(",", gpl@dataTable@table[, ind]))) {
			qqcat("multiple entries in one field\n")
		}
		df = data.frame("ID" = gpl@dataTable@table[, "ID"], "ENTREZID" = gpl@dataTable@table[, ind], stringsAsFactors = FALSE)

	} else if("GB_ACC" %in% all_nm) {
		qqcat("@{platform}: ID -> ACCNUM -> ENTREZID\n")
		ind = which(all_nm == "GB_ACC")
		if(any(grepl(",", gpl@dataTable@table[, ind]))) {
			qqcat("multiple entries in one field\n")
		}
		df = data.frame("ID" = gpl@dataTable@table[, "ID"], 
			"ACCNUM" = gsub("\\.\\d+$", "", gpl@dataTable@table[, ind]), 
			stringsAsFactors = FALSE)
		df$ENTREZID = ACCNUM2ENTREZID[df$ACCNUM]
		
	} else if("UNIGENE" %in% all_nm) {
		qqcat("@{platform}: ID -> UNIGENE -> ENTREZID\n")
		ind = which(all_nm == "UNIGENE")
		if(any(grepl(",", gpl@dataTable@table[, ind]))) {
			qqcat("multiple entries in one field\n")
		}
		df = data.frame("ID" = gpl@dataTable@table[, "ID"], "UNIGENE" = gpl@dataTable@table[, ind], stringsAsFactors = FALSE)
		if(!any(grepl("Hs", df$UNIGENE))) {
			df$UNIGENE = paste0("Hs.", df$UNIGENE)
		}
		df$ENTREZID = UNIGENE2ENTREZID[as.character(df$UNIGENE)]

	} else if(platform %in% c("GPL2986") && "GENE" %in% all_nm) {
		qqcat("@{platform}: ID -> ENTREZID\n")
		ind = which(all_nm == "GENE")
		if(any(grepl(",", gpl@dataTable@table[, ind]))) {
			qqcat("multiple entries in one field\n")
		}

		df = data.frame("ID" = gpl@dataTable@table[, "ID"], 
			"ENTREZID" = sapply(strsplit(gpl@dataTable@table[, ind], ","), function(x) {
				if(length(x) >= 1) {
					sample(x, 1)
				} else {
					NA
				}
			}), 
			stringsAsFactors = FALSE)

	} else if(platform %in% c("GPL3558") && "genbank" %in% all_nm) {
		qqcat("@{platform}: ID -> ACCNUM -> ENTREZID\n")
		ind = which(all_nm == "genbank")
		if(any(grepl(",", gpl@dataTable@table[, ind]))) {
			qqcat("multiple entries in one field\n")
		}
		df = data.frame("ID" = gpl@dataTable@table[, "ID"], "ACCNUM" = gpl@dataTable@table[, ind], stringsAsFactors = FALSE)
		df$ENTREZID = ACCNUM2ENTREZID[df$ACCNUM]

	} else if(platform %in% c("GPL2013") && "GENE" %in% all_nm) {
		qqcat("@{platform}: ID -> SYMBOL -> ENTREZID\n")
		ind = which(all_nm == "GENE")
		if(any(grepl(",", gpl@dataTable@table[, ind]))) {
			qqcat("multiple entries in one field\n")
		}
		df = data.frame("ID" = gpl@dataTable@table[, "ID"], "SYMBOL" = gpl@dataTable@table[, ind], stringsAsFactors = FALSE)
		df$ENTREZID = SYMBOL2ENTREZID[df$SYMBOL]

	} else if(platform %in% c("GPL3877") && "Gene Symbol" %in% all_nm) {
		qqcat("@{platform}: ID -> SYMBOL -> ENTREZID\n")
		ind = which(all_nm == "Gene Symbol")
		if(any(grepl(",", gpl@dataTable@table[, ind]))) {
			qqcat("multiple entries in one field\n")
		}
		df = data.frame("ID" = gpl@dataTable@table[, "ID"], "SYMBOL" = gpl@dataTable@table[, ind], stringsAsFactors = FALSE)
		df$ENTREZID = SYMBOL2ENTREZID[df$SYMBOL]
	} else {  ## GB_LIST
		qqcat("@{platform}: ID -> GB_LIST -> ENTREZID\n")
		ind = which(all_nm == "GB_LIST")
		if(any(grepl(",", gpl@dataTable@table[, ind]))) {
			qqcat("multiple entries in one field\n")
		}
		df = data.frame("ID" = gpl@dataTable@table[, "ID"], "GB_LIST" = gpl@dataTable@table[, ind], stringsAsFactors = FALSE)
		gb_list = strsplit(df$GB_LIST, ",")
		mapping = c(ACCNUM2ENTREZID, REFSEQ2ENTREZID[setdiff(names(REFSEQ2ENTREZID), names(ACCNUM2ENTREZID))])
		# mapping = hash(mapping)
		mapping = hashmap(names(mapping), mapping)
		counter = set_counter(length(gb_list), fmt = "converting GB_list to ENTREZID, %s")
		df$ENTREZID = sapply(gb_list, function(x) {
			counter()
			if(length(x) == 0) {
				return(NA)
			}
			x = gsub("\\.\\d+$", "", x)
			# y = sapply(x, function(x2) {
			# 	if(exists(x2, envir = mapping)) {
			# 		get(x2, envir = mapping)
			# 	} else {
			# 		NA
			# 	}
			# })
			y = mapping[[x]]
			y = unique(y)
			y = y[!is.na(y)]
			if(length(y > 1)) {
				return(sample(y, 1))
			} else if(length(y == 1)) {
				return(y)
			} else {
				return(NA)
			}
		})
	}

	df = df[, c("ID", "ENTREZID")]

	df = df[sample(nrow(df), nrow(df)), ]
	df = df[!duplicated(df[, 1]), ]
	
	map = structure(as.character(df[, 2]), names = df[, 1])
	qqcat("@{sum(!is.na(map))}/@{length(map)} have mappings\n")
	return(map)
}

all_platform = sapply(dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/", pattern = "^GPL\\d+\\.rds$"), function(id) gsub("\\.rds$", "", id))
for(platform in all_platform) {
	cat("==========================================\n")
	x = probe_id_2_entrez_id(platform)
	saveRDS(x, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/@{platform}_probe_id_to_entrez_id.rds"))
}

#########################################################
### for each of GDS datasets, try whether id mapping works
for(id in dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/", pattern = "GDS")) {
	# res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_all.rds"))
	gds = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}.rds"))
	
	# probe_id = rownames(get_matrix(res_list))
	probe_id = gds@dataTable@table$ID_REF
	platform = gds@header$platform
	id_mapping = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS_GPL/@{platform}_probe_id_to_entrez_id.rds"))

	entrez_id = id_mapping[probe_id]
	r = sum(!is.na(entrez_id))/length(entrez_id)
	if(r < 0.5)	qqcat("@{id}/@{platform}: mapped/all: @{sum(!is.na(entrez_id))}/@{length(entrez_id)}\n")
}

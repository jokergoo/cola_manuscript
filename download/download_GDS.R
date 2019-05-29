
# MESH terms for searching
# "Homo sapiens"[porgn] AND 50[n_samples] : 500[n_samples] AND "gds"[Filter]

id_list = scan("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/script/download/GDS_id_list.txt", what = "character")
library(GEOquery)
library(GetoptLong)

for(id in id_list) {
	if("GDS1761" %in% id) next
	dir.create(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}"))
	try({
		gds = getGEO(id, destdir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}"))
		saveRDS(gds, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}.rds"))
	})
}

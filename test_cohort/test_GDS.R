library(GetoptLong)

cores = 4
GetoptLong(
	"pid=s", "pid",
	"cores=i", "cores"
)

library(GEOquery)


######################################################
## get the expression matrix
######################################################
gds = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}.rds"))

ind = which(is.na(gds@dataTable@table[, 1]))
if(length(ind)) {
	qqcat("find @{length(ind)} rows have NA IDs.\n")
	gds@dataTable@table = gds@dataTable@table[-ind, ]
}
oe = try(eset <- GDS2eSet(gds))
if(inherits(oe, "try-error")) eset = GDS2eSet(gds, getGPL=FALSE)
mat = exprs(eset)

######################################################
## test whether need to perform log2 transformation
######################################################

## if in `value_type` column contain "log", most probably it means it is a two-channel
## microarray and the values are log fold change
if(!grepl("log", gds@header$value_type)) {
	mat[mat < 0] = 0
	x = mat[mat >= 0]
	ks_stat = ks.test(x, "pnorm", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))$statistic

	if(ks_stat > 0.3) {
		mat = log2(mat + 1)
	}
}

# samples with too many NA values
if(id == "GDS2808") {
	mat = mat[, apply(mat, 2, function(x) sum(is.na(x))/length(x) < 0.25), drop = FALSE]
}

library(cola)
register_NMF()
mat = adjust_matrix(mat)


################################################
## apply quantile normalization 
################################################
if(!grepl("log", gds@header$value_type)) {
	library(preprocessCore)
	cn = colnames(mat)
	rn = rownames(mat)
	mat = normalize.quantiles(mat)
	colnames(mat) = cn
	rownames(mat) = rn
}

##############################################
## get annotations
##############################################

anno = gds@dataTable@columns
anno = anno[, sapply(anno, function(x) {
	s = length(table(x))
	s > 1 && s < length(x)
}), drop = FALSE]
if(ncol(anno) == 0) anno = NULL

##############################################
## apply cola analysis
##############################################

res_list = run_all_consensus_partition_methods(mat, mc.cores = cores, anno = anno)
saveRDS(res_list, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_all.rds"))
cola_report(res_list, output_dir = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS/@{id}/@{id}_cola_report"), 
	mc.cores = cores, title = qq("cola Report for @{id}"))


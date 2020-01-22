
library(GetoptLong)
library(ComplexHeatmap)
library(circlize)

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test"


pmat1 = NULL
pmat2 = NULL

set.seed(123)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 5, nc = 6)))
i = 0
for(fraction in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
	i = i + 1
	j = 0
	membership_lt = lapply(1:100, function(x) c(rep(1, round(100*fraction)), rep(2, 100 - round(100*fraction))))

	p1 = NULL
	p2 = NULL
	for(p in c(seq(0, 0.1, by = 0.02))) {
		j = j + 1
		membership_lt2 = lapply(membership_lt, function(x) {
			l = sample(c(TRUE, FALSE), length(x), p = c(p, 1-p), replace = TRUE)
			x[l] = ifelse(x[l] == 1, 2, 1)
			x
		})
		membership_each = do.call("cbind", membership_lt2)

		consensus_mat = cola:::get_consensus_matrix(membership_each)

		pushViewport(viewport(layout.pos.row = i, layout.pos.col = j))
		pushViewport(viewport(width = unit(1, "npc") - unit(4, "mm"), height = unit(1, "npc") - unit(4, "mm")))
		# grid.rect()
		# f = ecdf(consensus_mat[lower.tri(consensus_mat)])
		# x = seq(0, 1, length = 100)
		# y = f(x)
		# x = c(0, x)
		# y = c(0, y)
		# grid.lines(x, y)
		# grid.lines(c(0.1, 0.1), c(0, 1), gp = gpar(lty = 2, col = "grey"))
		# grid.lines(c(0.9, 0.9), c(0, 1), gp = gpar(lty = 2, col = "grey"))
		
		p1 = c(p1, PAC_origin(consensus_mat))
		p2 = c(p2, aPAC(consensus_mat))

		ht = Heatmap(consensus_mat, col = colorRamp2(c(0, 1), c("white", "blue")),
			show_row_dend = FALSE, show_column_dend = FALSE, show_heatmap_legend = FALSE)
		draw(ht, newpage = FALSE, padding = unit(c(0, 0, 0, 0), "mm"))

		grid.text(qq("PAC_origin = @{sprintf('%.2f', p1[length(p1)])}\naPAC = @{sprintf('%.2f', p2[length(p2)])}"), 
			x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontsize = 6))

		popViewport()
		popViewport()
	}
	names(p1) = paste0("p_", c(seq(0, 0.1, by = 0.02)))
	names(p2) = paste0("p_", c(seq(0, 0.1, by = 0.02)))

	pmat1 = cbind(pmat1, p1)
	pmat2 = cbind(pmat2, p2)
}

colnames(pmat1) = paste0("frac_", c(0.1, 0.2, 0.3, 0.4, 0.5))
colnames(pmat2) = paste0("frac_", c(0.1, 0.2, 0.3, 0.4, 0.5))


plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "PAC_origin", ylab = "PAC_cola")
for(i in 1:ncol(pmat1)) {
	lines(pmat2[, i], pmat1[, i], col = i, type = "b")
}
abline(h = 0.1, lty = 2, col = "grey")
abline(v = 0.1, lty = 2, col = "grey")


##############################

project = "GDS"

collect_cm = function(project) {
	if(project == "recount2") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2")
	} else if(project == "GDS") {
		id_list = dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GDS", pattern = "^GDS\\d+$")
	}
	
	lt = list()
	for(id in id_list) {
		qqcat("loading @{id}...\n")
		oe = try({
			res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_cola_all.rds"))
			
		})

		if(!inherits(oe, "try-error")) {
			for(tm in res_list@top_value_method) {
				for(pm in res_list@partition_method) {
			
					res = res_list[tm, pm]
					if(!is.na(suggest_best_k(res))) {
						cm = get_consensus(res, k = suggest_best_k(res))
						attr(cm, "stat_foo") = stat_foo(res)
						lt[[paste(id, tm, pm, sep = "-")]] = cm
					}
				}
			}
		}
	}
	return(lt)
}

lt = collect_cm(project)

pac_cola = sapply(lt, cola:::PAC)
pac_origin = sapply(lt, PAC_origin)
pac_adapted = sapply(lt, aPAC)
pac_adapted_0.9 = sapply(lt, aPAC_foo, 0.9)
pac_adapted_0.8 = sapply(lt, aPAC_foo, 0.8)
pac_adapted_0.7 = sapply(lt, aPAC_foo, 0.7)
pac_foo = sapply(lt, function(x) attr(x, "stat_foo"))
pac_df = cbind(cola = pac_cola, origin = pac_origin, adapted = pac_adapted)
pac_df = pac_df[order(pac_df[, 2] - pac_df[, 1], decreasing = TRUE), ]

pdf(qq("@{BASE_DIR}/image/pac_compare.pdf"), width = 8, height = 8)
# plot(pac_df[, 2], pac_df[, 1], xlim = c(0, 1), ylim = c(0, 1), pch = 16, cex = 0.5, col = "#00000040")
# plot(pac_df[, 2], pac_df[, 3], xlim = c(0, 1), ylim = c(0, 1), pch = 16, cex = 0.5, col = "#00000040")
# plot(pac_df[, 2], pac_df[, 4], xlim = c(0, 1), ylim = c(0, 1), pch = 16, cex = 0.5, col = "#00000040")
plot(pac_adapted_0.9, pac_adapted, xlim = c(0, 1), ylim = c(0, 1), pch = 16, cex = 0.5, col = "#00000040")
plot(pac_adapted_0.8, pac_adapted, xlim = c(0, 1), ylim = c(0, 1), pch = 16, cex = 0.5, col = "#00000040")
plot(pac_adapted_0.7, pac_adapted, xlim = c(0, 1), ylim = c(0, 1), pch = 16, cex = 0.5, col = "#00000040")
dev.off()

pdf(qq("@{BASE_DIR}/image/pac_compare_heatmap.pdf"), width = 8, height = 8)
pac_df2 = pac_df[pac_df[, 1] < 0.1, ]
for(i in 1:20) {
	nm = rownames(pac_df2)[i]
	id = gsub("-.*$", "", nm)
	tm = gsub("^.*?-(.*?)-.*$", "\\1", nm)
	pm = gsub("^.*?-.*?-", "", nm)

	qqcat("loading @{id}...\n")
	res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/@{project}/@{id}/@{id}_cola_all.rds"))
	res = res_list[tm, pm]
	k = suggest_best_k(res)
	consensus_heatmap(res, k = k)
	membership_heatmap(res, k = k)
	plot_ecdf(res)
}
dev.off()


pdf(qq("@{BASE_DIR}/image/pac_compare_test.pdf"), width = 8, height = 8)
for(k in 2:6) {
	cm = get_consensus(res, k = k)
	plot(density(cm[lower.tri(cm)]))
}
plot_ecdf(res)
dev.off()





The scripts for analyzing GDS and recount2 datasets:

- `test_GDS.R`: R script that runs _cola_ analysis for every GDS dataset.
- `test_recount.R`: R script that runs _cola_ analysis for every recount dataset.
- `run_all_datasets.R`: R script that submits _cola_ analysis to computing clusters as well as hosts _cola_ reports on GitHub.

Downstream analysis:

- `compare_methods.R`: R script that makes plots in Figure 4 as well as in supplementaries.
- `enrichment.R`: R script that applies functional analysis.
- `id_mapping.R`: Because in GDS datasets, probe IDs are used as gene IDs (or row names of the matrix),
  to perform functional enrichment analysis, probe IDs need to be converted to Entrez IDs. This script
  tries to look for such information in the GPL annotation for each GDS dataset.


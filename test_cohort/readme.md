

- `id_mapping.R`: Because in GDS datasets, probe IDs are used as gene IDs (or row names of the matrix),
  to perform functional enrichment analysis, probe IDs need to be converted to Entrez IDs. This script
  tries to look for such information in the GPL annotations.
- `test_GDS.R`: run cola analysis for every GDS dataset.
- `test_recount.R`: run cola analysis for every recount dataset.
- `compare_methods.R`: 

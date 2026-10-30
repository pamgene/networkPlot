# networkPlot

Pathway-enrich and visualise a kinase network built by
[networkGen](https://github.com/pamgene/networkGen) (PCSF), or pathway-enrich
UKA data directly.

Two halves of one workflow:

- **`enrich_network()`** — pathway enrichment of the network's input kinases
  (`network_roles()`: connector nodes such as PCSF `Hidden` nodes or SIGNOR
  bridges are never tested), returning the pathway tables plus a per-node
  `pathway` annotation. **`enrich_uka()`** runs the same enrichment on UKA
  data directly, per comparison, on the input kinases networkGen's
  `uka_top()` selects. Both use **`enrich_kinases()`**:
  - pathway set: Reactome v2 (`reactome_v2_pathway_table()`, bundled in
    `inst/extdata`, see `data-raw/pathway-data.R`), with mutant/defect
    pathways that duplicate their normal counterpart left out
    (`filter_disease_pathways()`, judged with the run's kinase universe);
  - background: the UKAdb kinases (`ukadb_kinases()`), optionally restricted
    to the network's interactome, in at least one tested pathway;
  - test: one-sided Fisher + Benjamini-Hochberg; two tables:
    `pathways_<label>_all.csv` (every tested pathway with p < 0.1) and
    `pathways_<label>.csv` (the selected pathways, one per hit list: the
    smallest p, the others in `also_matches`).
- **`plot_network()`** — the interactive `visNetwork` HTML, with every style
  choice on **`networkplot_theme()`** (node shapes per type, the diverging
  `LogFC` fill, high-degree highlighting, and edge boldness mapped from
  interaction `cost` or run-frequency `weight`). `plot_pathway_heatmaps()`
  draws the kinase-by-pathway heatmaps: by default the "wrapped" layout (one
  small panel per pathway; one figure per comparison and one per pair of
  comparisons with their shared pathways; `pair_groups` limits the pairs to
  one group, e.g. a UKA Supergroup), or `layout = "classic"`.
  `plot_uka_pathway_heatmaps()` draws them from `enrich_uka()`'s result.

## Development

This project uses **R 4.3.0** (`/c/Program Files/R/R-4.3.0/bin/Rscript.exe`).
`ComplexHeatmap` and `circlize` (heatmaps) and `networkGen` (`enrich_uka()`)
are `Suggests` — the functions check for them at call time.

```r
devtools::test()
devtools::check()
```

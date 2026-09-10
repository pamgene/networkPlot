# networkPlot

Pathway-enrich and visualise a PCSF network built by
[networkGen](https://github.com/pamgene/networkGen).

Two halves of one workflow:

- **`enrich_network()`** — one Enrichr query (via the `enrichR` package)
  over every network node except the Steiner connectors (`type == "Hidden"`),
  returning a pathway table plus a per-node `pathway` annotation. Only
  Reactome is wired end to end today; the reference database is a pluggable
  `postprocess` step, and `download_reactome_refs()` refreshes the bundled
  Reactome snapshot.
- **`plot_network()`** — the interactive `visNetwork` HTML, with every style
  choice on **`networkplot_theme()`** (node shapes per type, the diverging
  `LogFC` fill, high-degree highlighting, and edge boldness mapped from
  interaction `cost` or run-frequency `weight`). `plot_pathway_heatmaps()`
  draws the kinase-by-pathway heatmaps.

`reconcile_pathways()` re-picks each comparison's pathways so the same
biology lines up across comparisons in the combined heatmaps.

## Development

This project uses **R 4.3.0** (`/c/Program Files/R/R-4.3.0/bin/Rscript.exe`).
`enrichR`, `ComplexHeatmap`, and `circlize` are `Suggests` — the enrichment
and heatmap functions check for them at call time.

```r
devtools::test()
devtools::check()
```

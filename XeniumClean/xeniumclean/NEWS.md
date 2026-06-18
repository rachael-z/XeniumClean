# XeniumClean 0.1.0

Initial release.

## Features

- `BuildReferenceGeneSets()`: build per-cell-type expressed/not-expressed gene
  sets from a single-cell RNA-seq reference, with configurable thresholds and
  ubiquitous-gene removal.
- `XeniumClean()`: main cleanup function. Works on a single Seurat object or a
  named list of Seurat objects (one per tissue section). Handles Seurat v4 and
  v5 assays, joins split layers transparently.
- `GetSpatialCoords()` / `ComputeSpatialNeighbors()`: flexible coordinate
  accessors and spatial neighbour computation.
- `CollapseLabels()`: helper for mapping fine spatial cluster labels to broader
  reference categories.
- `ScoreCellTypes()` / `FlagMismatches()`: per-cell scoring and identification
  of potentially misclassified cells.
- `PlotRemovedTranscripts()` / `CleanupSummary()`: QC visualisation and
  per-sample summaries.

## Known limitations

- Cleanup is conservative for cell types absent from the reference: those cells
  pass through untouched. Add custom references if you need them cleaned.
- Lineage-continuum cell types (T/NK, B/Plasma, Mono/Mac) may benefit from the
  `skip.pairs` argument to avoid over-cleaning between closely related types.

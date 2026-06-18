# XeniumClean

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Spatial neighbour-aware transcript cleanup for imaging-based spatial transcriptomics data (Xenium, MERSCOPE, CosMx).

XeniumClean removes biologically implausible transcripts by combining a single-cell RNA-seq reference with spatial neighbourhood information. For every cell, it identifies which neighbouring cells belong to other cell types and erases only those transcripts that the reference says cannot plausibly originate from the cell's own type.

## How it works

For each cell type, two gene sets are built from a single-cell reference:

- **Expressed genes**: genes detected in at least 10% of cells of that type
- **Not-expressed genes**: genes detected in at most 5% of cells

Then, for each cell `x` of type A in spatial data:

1. Find spatial neighbours within a chosen radius (default 50 µm)
2. For each neighbouring cell type B (B != A), genes that are *expressed in B but not in A* are biologically implausible in A
3. If those genes appear in cell `x`'s transcript counts, erase them
4. Never remove genes that A could genuinely express

The output is a new cleaned assay (`XeniumClean` by default) alongside the original counts.

## ⚠️ Critical warnings — read before running

### 1. Cell types must be pre-defined in both datasets

Before running XeniumClean, both the single-cell reference and your spatial object must already have a cell type column populated. XeniumClean does not perform clustering or annotation — it only acts on pre-existing labels.

If the cell type names differ between reference and spatial (e.g. reference has `"T cells CD4+"`, your spatial data has `"CD4_T_cells"`), use `CollapseLabels()` to align them before calling `XeniumClean()`. The names must match exactly for the gene set lookup to work.

### 2. Use BROAD cell types, not fine subsets

Define cell types at the **major lineage level** (e.g. `T_cells`, `B_cells`, `Macrophage`, `Cancer`), NOT at the subset level (e.g. `CD4_T_cells`, `Memory_B_cells`, `M2_Macrophage`).

Why: closely related subsets share most of their expression program. If you split T cells into CD4 vs CD8 for cleanup, genes that are real in both (TCF7, IL7R, RUNX3) may end up in CD4's `NotExpressedGenes` for the CD8 reference and vice versa, leading to over-cleaning. Even grouping NK cells with T cells (as "lymphoid") is reasonable if your panel can't reliably tell them apart.

**The rule**: when in doubt, group more coarsely.

### 3. Post-cleanup is *better* for subclustering

This is actually a major benefit of XeniumClean. In raw Xenium data, fine subclustering often picks up neighbourhood patterns rather than true biology — T cells near macrophages cluster together because they're contaminated with macrophage genes, not because they're a "macrophage-adjacent T cell" subset. After cleanup, those neighbourhood artifacts are gone, so subclustering on the cleaned assay reveals genuine subsets (CD4 vs CD8, memory vs effector, etc.) based on real expression.

**Recommended workflow:**
- Cluster broadly → label cell types (BROAD: T_cells, B_cells, etc.)
- Run XeniumClean using these broad labels
- THEN subcluster on the cleaned assay to find true subsets

### 4. ⚠️ Reference choice is critical — wrong references erase biology

XeniumClean is only as good as its reference. A mismatched reference can:

- Mark genes as "not expressed" that ARE expressed in your tissue → those genes get erased = real biology destroyed
- Mark genes as "expressed" when they shouldn't be → contaminating transcripts kept = no cleanup happens

**Recommendations, in order of safety:**

1. ✅ **Best**: single-cell data from the same experiment (e.g. matched scRNA-seq from a subset of your samples)
2. ✅ **Good**: a published atlas from the same tissue, disease state, and sample preparation (e.g. Wu 2021 breast cancer for breast cancer Xenium)
3. ⚠️ **Risky**: a public atlas from a different tissue or different disease (e.g. healthy breast for breast cancer; PBMC for tumor)
4. ❌ **Avoid**: manually curated gene lists from literature alone

For option 1 or 2, the gene sets will broadly reflect your tissue's biology. For option 3, expect over-cleaning of disease-specific or tissue-specific genes that the reference doesn't capture.

**User-defined references**: you can manually construct `gene.sets` with `ExpressedGenes` and `NotExpressedGenes` lists, but this is extremely risky. You'd need to accurately know, for every cell type and every gene in your panel, whether it's expressed (≥10%) or not (≤5%). For a 380-gene panel × 10 cell types, that's 3800 gene-level decisions — almost impossible to get right manually. Use a single-cell reference unless you have a very specific reason not to.

### 5. Multi-sample Seurat objects must be split into a list

XeniumClean operates per-section because neighbours are within-section. If you have a merged Seurat object with multiple samples, you must split it into a named list before passing to `XeniumClean()`.

**Do not use `subset()`** for v5 objects with split layers — it can be slow and produce inconsistent results. Use one of:

```r
# Option 1: SplitObject (Seurat built-in)
xenium_list <- SplitObject(xenium_merged, split.by = "sample")

# Option 2: subset_opt (community function, faster for v5 split-layer objects)
# See https://github.com/satijalab/seurat/issues/2724
xenium_list <- lapply(unique(xenium_merged$sample), function(s) {
  subset_opt(xenium_merged, subset = sample == s)
})
names(xenium_list) <- unique(xenium_merged$sample)
```

The split must produce one Seurat object per section, each with its own FOV (or spatial coordinates), so that neighbours can be computed within section.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rzemek/XeniumClean")
```

## Quick start

```r
library(XeniumClean)
library(Seurat)

# 1. Build per-cell-type gene sets from a single-cell reference
shared_genes <- intersect(rownames(sc_reference), rownames(xenium_obj))

gene_sets <- BuildReferenceGeneSets(
  reference = sc_reference,
  group.by  = "cell_type_broad",          # BROAD labels (T_cells, not CD8_T_cells)
  genes.use = shared_genes,
  expressed.threshold     = 0.10,
  not.expressed.threshold = 0.05
)

# 2. (Optional) Collapse fine spatial labels to match broad reference granularity
xenium_obj$xclean_label <- CollapseLabels(
  labels  = xenium_obj$clusters,
  mapping = c(
    CD4_T_cells   = "T_cells",
    CD8_T_cells   = "T_cells",
    T_reg_cells   = "T_cells",
    NK_cells      = "T_cells",            # consider grouping with T_cells
    Cancer_LumA   = "Cancer",
    Cancer_Basal  = "Cancer"
  )
)

# 3. Validate setup before running (recommended)
ValidateSetup(xenium_obj, gene_sets, group.by = "xclean_label")

# 4. Run cleanup
xenium_obj <- XeniumClean(
  object         = xenium_obj,
  gene.sets      = gene_sets,
  group.by       = "xclean_label",
  radius         = 50,
  new.assay.name = "XeniumClean"
)

# 5. QC
CleanupSummary(xenium_obj, group.by = "xclean_label")
PlotRemovedTranscripts(xenium_obj, group.by = "xclean_label")

# 6. Use the cleaned assay for downstream analysis - INCLUDING SUBCLUSTERING
DefaultAssay(xenium_obj) <- "XeniumClean"
xenium_obj <- NormalizeData(xenium_obj)
xenium_obj <- ScaleData(xenium_obj)
xenium_obj <- RunPCA(xenium_obj, features = rownames(xenium_obj))
xenium_obj <- RunUMAP(xenium_obj, dims = 1:30)

# Subclustering on the cleaned assay reveals true subsets, not contamination patterns
xenium_obj <- FindNeighbors(xenium_obj, dims = 1:30)
xenium_obj <- FindClusters(xenium_obj, resolution = 0.8)
```

## Working with multiple sections

```r
# Split a merged object into per-sample list first
xenium_list <- SplitObject(xenium_merged, split.by = "sample")

# Then run cleanup across the list
xenium_list <- XeniumClean(
  object    = xenium_list,
  gene.sets = gene_sets,
  group.by  = "xclean_label",
  radius    = 50,
  parallel  = TRUE,
  workers   = 6
)
```

## Different object layouts

XeniumClean tries to find coordinates automatically (image FOVs, then `spatial` reduction, then metadata columns). If you need to override, pass the source:

```r
# Coordinates in a 'spatial' dimensional reduction
XeniumClean(xenium_obj, ..., coords.source = "reduction", reduction.name = "spatial")

# Coordinates in custom metadata columns
XeniumClean(xenium_obj, ..., coords.source = "metadata", coord.cols = c("cell_x","cell_y"))

# Or pass coordinates explicitly
XeniumClean(xenium_obj, ..., coords.source = "matrix", coords = my_coord_matrix)
```

`XeniumClean` also works with both Seurat v4 (`slot = "counts"`) and v5 (`layer = "counts"`) objects, including v5 objects with split layers (it joins them transiently without modifying the original).

## Tuning

| Parameter | Default | When to change |
|---|---|---|
| `expressed.threshold` | 0.10 | Lower to be more aggressive (more genes considered "expressed"); raise for stricter cell-type-specificity |
| `not.expressed.threshold` | 0.05 | Lower to be more conservative (fewer removals); raise to remove more |
| `radius` | 50 µm | Distance between centroids of neighbouring cells. Match your imaging resolution |
| `skip.pairs` | `NULL` | Add closely-related lineages (T/NK, B/Plasma, Mono/Mac) to prevent over-cleaning |
| `remove.ubiquitous` | `TRUE` | Set `FALSE` to keep housekeeping genes in cell-type marker sets |
| `block.size` | 500 | Reduce for memory-constrained workers, increase for fewer overhead calls |

## Documentation

- **Vignette**: `vignette("xeniumclean")` — full walkthrough with a worked example
- **Function reference**: `?XeniumClean`, `?BuildReferenceGeneSets`, etc.

## Citation

If you use XeniumClean in your research, please cite:

> Zemek RM. XeniumClean: spatial neighbour-aware transcript cleanup for imaging-based
> spatial transcriptomics. https://github.com/rzemek/XeniumClean

## License
OPEN

#' XeniumClean: Spatial neighbour-aware transcript cleanup
#'
#' XeniumClean removes biologically implausible transcripts from imaging-based
#' spatial transcriptomics data (Xenium, MERSCOPE, CosMx) by combining a
#' single-cell RNA-seq reference with spatial neighbourhood information.
#'
#' @section Concept summary:
#' For each cell type, the package defines two gene sets from a single-cell
#' reference:
#' - **Expressed genes**: genes detected in a meaningful fraction of cells
#'   (default >=10%)
#' - **Not-expressed genes**: genes detected in a small fraction of cells
#'   (default <=5%)
#'
#' Then, for each cell `x` of type A in the spatial data:
#' 1. Identify its spatial neighbours within a chosen radius.
#' 2. For each neighbouring cell type B that differs from A, find genes that
#'    are expressed in B but not expressed in A.
#' 3. Those genes are biologically implausible in A, so if they appear in cell
#'    `x`'s transcript counts, set them to zero ("erase").
#' 4. Never remove genes that A could genuinely express.
#'
#' This produces spatial, biologically consistent decontamination that respects
#' physical adjacency rather than relying on global expression patterns alone.
#'
#' @section Typical workflow:
#' 1. [BuildReferenceGeneSets()] from a published scRNA-seq atlas
#' 2. [CollapseLabels()] to align reference and spatial cell type granularity
#' 3. [XeniumClean()] to perform the actual cleanup
#' 4. [PlotRemovedTranscripts()] / [CleanupSummary()] for QC
#' 5. (Optional) [ScoreCellTypes()] / [FlagMismatches()] for misclassification
#'    diagnosis
#'
#' @section Design philosophy:
#' - Works on Seurat v4 and v5 objects
#' - Handles coordinates stored in image FOVs, spatial reductions, or metadata
#' - Tunable thresholds, radius, and parallelism
#' - Conservative: only removes transcripts where there is direct spatial
#'   evidence of a different cell type as the source
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom methods is
NULL

# Avoid R CMD check NOTE about .data
utils::globalVariables(c(".data"))

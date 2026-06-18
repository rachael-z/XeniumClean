#' Build per-cell-type gene sets from a single-cell reference
#'
#' For each cell type in a reference Seurat object, computes:
#' - `ExpressedGenes`: genes expressed (non-zero) in at least
#'   `expressed.threshold` fraction of cells of that type.
#' - `NotExpressedGenes`: genes expressed in at most `not.expressed.threshold`
#'   fraction.
#'
#' Genes between the two thresholds are considered ambiguous and are not used
#' for cleanup. Genes expressed across all cell types ("ubiquitous") are
#' optionally removed from `ExpressedGenes` since they cannot discriminate
#' between cell types.
#'
#' The output is intended to be passed to [XeniumClean()] for the actual
#' transcript cleanup.
#'
#' @section Critical: reference choice:
#' \strong{XeniumClean is only as good as its reference.} A mismatched
#' reference can erase real biology (genes that ARE expressed in your tissue
#' get put in `NotExpressedGenes`) or fail to clean (genes that AREN'T
#' expressed get put in `ExpressedGenes`).
#'
#' In order of safety, the reference should be:
#' \enumerate{
#'   \item \strong{Best}: scRNA-seq from the same experiment (matched tissue,
#'     disease, sample prep)
#'   \item \strong{Good}: published atlas of the same tissue + disease
#'   \item \strong{Risky}: published atlas of healthy tissue when you have
#'     disease, or related-but-different tissue
#'   \item \strong{Avoid}: manually curated literature gene lists alone
#' }
#'
#' Always sanity-check the resulting gene sets against established markers:
#' `"CD3D" %in% gene_sets$ExpressedGenes$T_cells` should be `TRUE`, etc.
#'
#' @section Critical: cell type granularity:
#' Use \strong{broad} cell type categories (T_cells, B_cells, Macrophage,
#' Cancer), not fine subsets (CD4_T_cells, Memory_B_cells, M2_Macrophage).
#' Closely related subsets share most of their expression program; splitting
#' them puts shared markers in each other's `NotExpressedGenes`, leading to
#' over-cleaning. Subclustering for true subsets should be done on the cleaned
#' data, not before cleanup.
#'
#' @param reference A Seurat object containing the single-cell reference (e.g.
#'   a published breast cancer atlas, or a Flex/Chromium dataset).
#' @param group.by Metadata column on `reference` holding cell type labels.
#' @param assay Assay to pull normalised data from. Defaults to the reference's
#'   `DefaultAssay`.
#' @param layer Layer/slot holding log-normalised data (default `"data"`).
#' @param genes.use Optional character vector restricting the universe of genes
#'   considered. Use this to restrict to the Xenium panel:
#'   `intersect(rownames(reference), rownames(xenium_object))`.
#' @param expressed.threshold Fraction of cells (within a cell type) that must
#'   express a gene for it to be "expressed". Default `0.10`.
#' @param not.expressed.threshold Maximum fraction of cells that may express a
#'   gene for it to be "not expressed". Default `0.05`.
#' @param min.cells.per.group Minimum cells in a group to compute gene sets.
#'   Smaller groups are skipped with a warning. Default `30`.
#' @param remove.ubiquitous If `TRUE` (default), drop genes that appear in
#'   `ExpressedGenes` for every cell type.
#' @param verbose Logical; print progress per group.
#'
#' @return A list of class `XeniumCleanGeneSets` with elements:
#'   - `ExpressedGenes`: named list of character vectors (one per cell type)
#'   - `NotExpressedGenes`: named list of character vectors
#'   - `ubiquitous_genes`: genes dropped as ubiquitous
#'   - `genes_used`: the gene universe considered
#'   - `thresholds`: numeric vector with the two thresholds used
#'   - `cell_counts`: number of cells per group
#'
#' @examples
#' \dontrun{
#' # Restrict to genes shared with your Xenium panel
#' shared <- intersect(rownames(sc_reference), rownames(xenium_obj))
#'
#' gene_sets <- BuildReferenceGeneSets(
#'   reference  = sc_reference,
#'   group.by   = "celltype_minor",
#'   genes.use  = shared,
#'   expressed.threshold     = 0.10,
#'   not.expressed.threshold = 0.05
#' )
#'
#' # Sanity check
#' "CD3D"  %in% gene_sets$ExpressedGenes$T_cells       # TRUE
#' "CD3D"  %in% gene_sets$NotExpressedGenes$Cancer     # TRUE
#' }
#'
#' @export
BuildReferenceGeneSets <- function(reference,
                                   group.by,
                                   assay                    = NULL,
                                   layer                    = "data",
                                   genes.use                = NULL,
                                   expressed.threshold      = 0.10,
                                   not.expressed.threshold  = 0.05,
                                   min.cells.per.group      = 30,
                                   remove.ubiquitous        = TRUE,
                                   verbose                  = TRUE) {

  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(reference)
  .validate_object(reference, assay, group.by)

  if (expressed.threshold <= not.expressed.threshold) {
    stop("`expressed.threshold` must be greater than `not.expressed.threshold`.")
  }

  expr_mat <- .get_layer_data(reference, assay = assay, layer = layer)

  # Restrict to shared gene universe if provided
  if (!is.null(genes.use)) {
    genes.use <- intersect(genes.use, rownames(expr_mat))
    if (length(genes.use) == 0) {
      stop("No overlap between `genes.use` and reference rownames.")
    }
    expr_mat <- expr_mat[genes.use, ]
  } else {
    genes.use <- rownames(expr_mat)
  }

  labels <- as.character(reference@meta.data[[group.by]])
  celltypes <- sort(unique(labels[!is.na(labels)]))

  ExpressedGenes    <- list()
  NotExpressedGenes <- list()
  cell_counts       <- integer(0)

  for (cl in celltypes) {
    cells_in <- which(labels == cl)
    if (length(cells_in) < min.cells.per.group) {
      .v_message("Skipping ", cl, " (only ",
                 length(cells_in), " cells)", verbose = verbose)
      next
    }

    frac_in <- Matrix::rowSums(expr_mat[, cells_in, drop = FALSE] > 0) /
                length(cells_in)

    ExpressedGenes[[cl]]    <- names(frac_in[frac_in >= expressed.threshold])
    NotExpressedGenes[[cl]] <- names(frac_in[frac_in <= not.expressed.threshold])
    cell_counts[[cl]]       <- length(cells_in)

    .v_message(sprintf("%-25s n=%6d  expressed=%4d  not_expressed=%4d",
                       cl, length(cells_in),
                       length(ExpressedGenes[[cl]]),
                       length(NotExpressedGenes[[cl]])),
               verbose = verbose)
  }

  # Remove ubiquitous genes (in ExpressedGenes for every cell type)
  ubiquitous <- character(0)
  if (isTRUE(remove.ubiquitous) && length(ExpressedGenes) > 0) {
    counts <- table(unlist(ExpressedGenes))
    ubiquitous <- names(counts[counts == length(ExpressedGenes)])
    ExpressedGenes <- lapply(ExpressedGenes, function(g) setdiff(g, ubiquitous))
    .v_message("Removed ", length(ubiquitous),
               " ubiquitous genes from ExpressedGenes.",
               verbose = verbose)
  }

  structure(
    list(
      ExpressedGenes    = ExpressedGenes,
      NotExpressedGenes = NotExpressedGenes,
      ubiquitous_genes  = ubiquitous,
      genes_used        = genes.use,
      thresholds        = c(expressed     = expressed.threshold,
                            not_expressed = not.expressed.threshold),
      cell_counts       = cell_counts,
      reference_label   = group.by
    ),
    class = "XeniumCleanGeneSets"
  )
}


#' Print method for XeniumCleanGeneSets
#'
#' @param x A `XeniumCleanGeneSets` object.
#' @param ... Ignored.
#' @export
print.XeniumCleanGeneSets <- function(x, ...) {
  cat("XeniumCleanGeneSets\n")
  cat("  Cell types: ", length(x$ExpressedGenes), "\n", sep = "")
  cat("  Gene universe: ", length(x$genes_used), " genes\n", sep = "")
  cat("  Thresholds: expressed >= ", x$thresholds["expressed"],
      ", not expressed <= ", x$thresholds["not_expressed"], "\n", sep = "")
  cat("  Ubiquitous genes removed: ", length(x$ubiquitous_genes), "\n", sep = "")
  cat("  Per-cell-type counts:\n")
  for (cl in names(x$ExpressedGenes)) {
    cat(sprintf("    %-25s  expressed=%4d  not_expressed=%4d\n",
                cl,
                length(x$ExpressedGenes[[cl]]),
                length(x$NotExpressedGenes[[cl]])))
  }
  invisible(x)
}

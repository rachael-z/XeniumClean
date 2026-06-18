#' Validate setup before running XeniumClean
#'
#' Performs sanity checks before the actual cleanup. Catches the most common
#' setup errors:
#'
#' - Cell type label mismatch between spatial and reference
#' - Gene sets that fail basic biological sanity (e.g. CD3D not in T_cells)
#' - Coordinate availability
#' - Too-fine cell type granularity that XeniumClean is unlikely to handle well
#'
#' Run this once before your first `XeniumClean()` call. It returns invisibly
#' but prints a structured report.
#'
#' @param object A Seurat object or named list of Seurat objects.
#' @param gene.sets A `XeniumCleanGeneSets` object.
#' @param group.by Metadata column with cell type labels (same one you'd pass
#'   to `XeniumClean()`).
#' @param coords.source,fov,reduction.name,coord.cols Passed to
#'   [GetSpatialCoords()] for the coordinate check.
#' @param sanity.markers A named list of canonical markers per cell type, used
#'   to verify the gene sets are biologically plausible. Defaults to a small
#'   immune/stromal/cancer panel. Names must match cell type labels.
#'
#' @return Invisibly returns `TRUE` if all checks pass. Prints warnings and
#'   errors for any failures.
#'
#' @examples
#' \dontrun{
#' ValidateSetup(xenium_list, gene_sets, group.by = "xclean_label")
#' }
#'
#' @export
ValidateSetup <- function(object,
                          gene.sets,
                          group.by,
                          coords.source  = "auto",
                          fov            = NULL,
                          reduction.name = "spatial",
                          coord.cols     = c("x", "y"),
                          sanity.markers = NULL) {

  if (is.null(sanity.markers)) {
    sanity.markers <- list(
      T_cells      = c("CD3D","CD3E","CD2"),
      B_cells      = c("MS4A1","CD79A","CD19"),
      NK_cells     = c("NKG7","KLRD1","GNLY"),
      Plasma_cells = c("JCHAIN","MZB1","XBP1"),
      Macrophage   = c("CD68","CD163","C1QB"),
      Monocytes    = c("FCGR3A","S100A9","CD14"),
      DC           = c("CD1C","ITGAX","CD74"),
      Cancer       = c("EPCAM","KRT8","KRT18"),
      Fibroblasts  = c("DCN","LUM","COL1A1"),
      Endothelial  = c("PECAM1","CDH5","VWF")
    )
  }

  pass <- TRUE
  msg  <- function(level, txt) {
    cat(sprintf("  [%s] %s\n", level, txt))
    if (level == "FAIL") pass <<- FALSE
  }

  cat("=== XeniumClean setup check ===\n")
  .validate_gene_sets(gene.sets)

  # 1. group.by exists in the spatial object(s)
  cat("\n[1/5] group.by column present\n")
  if (methods::is(object, "Seurat")) objs <- list(.x = object) else objs <- object
  for (nm in names(objs)) {
    obj <- objs[[nm]]
    if (group.by %in% colnames(obj@meta.data)) {
      msg("OK", sprintf("'%s' present in %s", group.by, nm))
    } else {
      msg("FAIL", sprintf("'%s' NOT found in %s. Available: %s",
                          group.by, nm,
                          paste(colnames(obj@meta.data), collapse = ", ")))
    }
  }

  # 2. Label overlap between spatial and gene.sets
  cat("\n[2/5] Cell type label overlap with reference\n")
  spatial_labels <- unique(unlist(lapply(objs, function(o) {
    if (group.by %in% colnames(o@meta.data)) {
      as.character(o@meta.data[[group.by]])
    } else NA_character_
  })))
  spatial_labels <- spatial_labels[!is.na(spatial_labels)]
  ref_labels <- names(gene.sets$ExpressedGenes)

  matched   <- intersect(spatial_labels, ref_labels)
  unmatched <- setdiff(spatial_labels, ref_labels)
  ref_extra <- setdiff(ref_labels,    spatial_labels)

  msg("OK", sprintf("%d / %d spatial labels matched to reference",
                    length(matched), length(spatial_labels)))
  if (length(unmatched) > 0) {
    msg("WARN", sprintf("Spatial-only labels (will be SKIPPED, not cleaned): %s",
                        paste(unmatched, collapse = ", ")))
  }
  if (length(ref_extra) > 0) {
    msg("INFO", sprintf("Reference-only labels (unused): %s",
                        paste(ref_extra, collapse = ", ")))
  }
  if (length(matched) == 0) {
    msg("FAIL", "No label overlap. Use CollapseLabels() to align names.")
  }

  # 3. Granularity check
  cat("\n[3/5] Cell type granularity\n")
  finegrain_hints <- c("CD4", "CD8", "Treg", "Cycling", "Naive", "Memory",
                       "M1", "M2", "cDC1", "cDC2", "pDC", "LumA", "LumB",
                       "Basal", "Her2")
  finegrain_hits  <- vapply(matched, function(lbl) {
    any(vapply(finegrain_hints, function(h) grepl(h, lbl, ignore.case = TRUE),
               logical(1)))
  }, logical(1))
  if (any(finegrain_hits)) {
    msg("WARN",
        sprintf("Possible fine-subset labels detected: %s. Use BROAD labels for cleanup; subcluster AFTER cleanup.",
                paste(matched[finegrain_hits], collapse = ", ")))
  } else {
    msg("OK", "All matched labels look broad-level")
  }

  # 4. Biological sanity check on gene sets
  cat("\n[4/5] Biological sanity check on gene sets\n")
  for (cl in intersect(names(sanity.markers), names(gene.sets$ExpressedGenes))) {
    markers <- sanity.markers[[cl]]
    markers <- intersect(markers, gene.sets$genes_used %||% character(0))
    if (length(markers) == 0) next

    present <- markers %in% gene.sets$ExpressedGenes[[cl]]
    if (any(present)) {
      msg("OK", sprintf("%s: %s in ExpressedGenes",
                        cl, paste(markers[present], collapse = ", ")))
    } else {
      msg("WARN",
          sprintf("%s: NONE of %s found in ExpressedGenes (expected at least one). Check label mapping or reference quality.",
                  cl, paste(markers, collapse = ", ")))
    }
  }

  # 5. Coordinate availability
  cat("\n[5/5] Spatial coordinate availability\n")
  for (nm in names(objs)) {
    obj <- objs[[nm]]
    ok <- tryCatch({
      c <- GetSpatialCoords(obj,
                            coords.source  = coords.source,
                            fov            = fov,
                            reduction.name = reduction.name,
                            coord.cols     = coord.cols)
      nrow(c) >= 0.9 * ncol(obj)
    }, error = function(e) {
      msg("FAIL", sprintf("%s: %s", nm, conditionMessage(e)))
      FALSE
    })
    if (isTRUE(ok)) msg("OK", sprintf("%s: coordinates available", nm))
  }

  cat("\n=== ", if (pass) "ALL CHECKS PASSED" else "CHECKS FAILED — fix above before running XeniumClean()", " ===\n", sep = "")
  invisible(pass)
}

# null-coalescing operator (used above)
`%||%` <- function(a, b) if (is.null(a)) b else a

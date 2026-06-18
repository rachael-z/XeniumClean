#' Per-cell scoring against reference gene sets
#'
#' For each cell type in `gene.sets`, computes a score per cell defined as
#' `mean(positive markers) - mean(negative markers)` using log-normalised
#' expression. This penalises cells expressing wrong-cell-type markers, not
#' just rewards presence of correct ones.
#'
#' Adds `pred_top1`, `pred_top2`, `pred_score1`, `pred_score2`, and
#' `pred_margin` to the object metadata.
#'
#' @param object A Seurat object with a normalised data layer in `assay`.
#' @param gene.sets A `XeniumCleanGeneSets` from [BuildReferenceGeneSets()].
#' @param assay Source assay. Default `"Xenium"`.
#' @param layer Source layer (log-normalised data). Default `"data"`.
#' @param min.genes Minimum positive/negative markers required for a cell type
#'   to be scored (otherwise returns zero score for that type). Default `3`.
#'
#' @return The Seurat object with prediction columns added.
#'
#' @examples
#' \dontrun{
#' xenium_obj <- ScoreCellTypes(xenium_obj, gene_sets, assay = "Xenium")
#' head(xenium_obj@meta.data[, c("pred_top1","pred_top2","pred_margin")])
#' }
#'
#' @export
ScoreCellTypes <- function(object,
                           gene.sets,
                           assay     = "Xenium",
                           layer     = "data",
                           min.genes = 3) {

  .validate_gene_sets(gene.sets)
  norm_data <- .get_layer_data(object, assay = assay, layer = layer)

  score_mat <- vapply(
    names(gene.sets$ExpressedGenes),
    function(ct) {
      pos_genes <- intersect(gene.sets$ExpressedGenes[[ct]], rownames(norm_data))
      neg_genes <- intersect(gene.sets$NotExpressedGenes[[ct]], rownames(norm_data))

      pos_score <- if (length(pos_genes) >= min.genes) {
        Matrix::colMeans(norm_data[pos_genes, , drop = FALSE])
      } else rep(0, ncol(norm_data))

      neg_score <- if (length(neg_genes) >= min.genes) {
        Matrix::colMeans(norm_data[neg_genes, , drop = FALSE])
      } else rep(0, ncol(norm_data))

      pos_score - neg_score
    },
    FUN.VALUE = numeric(ncol(norm_data))
  )
  rownames(score_mat) <- colnames(norm_data)

  # Top-1 and top-2 prediction
  top1_idx   <- max.col(score_mat, ties.method = "first")
  top1_name  <- colnames(score_mat)[top1_idx]
  top1_score <- score_mat[cbind(seq_len(nrow(score_mat)), top1_idx)]

  score2 <- score_mat
  score2[cbind(seq_len(nrow(score2)), top1_idx)] <- -Inf
  top2_idx   <- max.col(score2, ties.method = "first")
  top2_name  <- colnames(score_mat)[top2_idx]
  top2_score <- score_mat[cbind(seq_len(nrow(score_mat)), top2_idx)]

  object$pred_top1   <- top1_name
  object$pred_top2   <- top2_name
  object$pred_score1 <- top1_score
  object$pred_score2 <- top2_score
  object$pred_margin <- top1_score - top2_score

  object
}


#' Flag potentially misclassified cells
#'
#' Identifies cells whose assigned label is neither the top-1 nor top-2 best
#' match by per-cell scoring. Optionally compares to spatial neighbour types
#' to distinguish contamination (neighbours match the predicted type, so the
#' cell is a real assigned-type cell with contamination) from misclassification
#' (neighbours don't match, suggesting the assignment itself is wrong).
#'
#' Adds `mismatch_flag` and optionally `xeniumclean_diagnosis` to metadata.
#'
#' @param object A Seurat object that has already been scored by
#'   [ScoreCellTypes()].
#' @param group.by Metadata column with assigned labels (must overlap with the
#'   gene sets used in scoring).
#' @param check.neighbors If `TRUE` (default), compute the spatial neighbour
#'   composition for each flagged cell and label the diagnosis as
#'   `"likely_contamination"`, `"ambiguous"`, or `"likely_misclassification"`.
#' @param radius Spatial radius for neighbour check (only used if
#'   `check.neighbors = TRUE`). Default `50`.
#' @param contamination.threshold A flagged cell is called
#'   `"likely_contamination"` if at least this fraction of its neighbours match
#'   the predicted type. Default `0.5`.
#' @param assigned.threshold A flagged cell is called `"ambiguous"` if at least
#'   this fraction of its neighbours match the assigned (original) type.
#'   Default `0.25`.
#' @param ... Passed to [GetSpatialCoords()] if neighbour check is enabled.
#'
#' @return The object with new metadata columns.
#'
#' @examples
#' \dontrun{
#' xenium_obj <- ScoreCellTypes(xenium_obj, gene_sets)
#' xenium_obj <- FlagMismatches(xenium_obj, group.by = "cell_type")
#' table(xenium_obj$xeniumclean_diagnosis, xenium_obj$cell_type)
#' }
#'
#' @export
FlagMismatches <- function(object,
                           group.by,
                           check.neighbors          = TRUE,
                           radius                   = 50,
                           contamination.threshold  = 0.5,
                           assigned.threshold       = 0.25,
                           ...) {

  if (!"pred_top1" %in% colnames(object@meta.data)) {
    stop("Run ScoreCellTypes() before FlagMismatches().")
  }
  if (!group.by %in% colnames(object@meta.data)) {
    stop("Metadata column '", group.by, "' not found.")
  }

  current      <- as.character(object@meta.data[[group.by]])
  predictable  <- unique(c(object$pred_top1, object$pred_top2))
  evaluable    <- current %in% predictable

  object$mismatch_flag <- evaluable &
    (current != object$pred_top1) &
    (current != object$pred_top2)

  if (!isTRUE(check.neighbors)) return(object)

  # Need neighbours
  if (is.null(object@misc$xeniumclean_neighbors)) {
    object <- ComputeSpatialNeighbors(object, radius = radius, ...)
  }
  neighbors_list <- object@misc$xeniumclean_neighbors

  diag_vec <- rep(NA_character_, ncol(object))
  flagged  <- which(object$mismatch_flag)
  if (length(flagged) == 0) {
    object$xeniumclean_diagnosis <- diag_vec
    return(object)
  }

  for (i in flagged) {
    nidx <- neighbors_list[[i]]
    if (length(nidx) == 0) next
    n_match_pred <- mean(current[nidx] == object$pred_top1[i])
    n_match_orig <- mean(current[nidx] == current[i])
    diag_vec[i] <- if (n_match_pred >= contamination.threshold &&
                       n_match_orig <  assigned.threshold) {
      "likely_contamination"
    } else if (n_match_orig >= assigned.threshold) {
      "ambiguous"
    } else {
      "likely_misclassification"
    }
  }

  object$xeniumclean_diagnosis <- diag_vec
  object
}


#' Apply a label collapse mapping
#'
#' Convenience wrapper around named-vector recoding. Maps fine cluster labels
#' to broader categories that match a reference's cell type granularity.
#'
#' @param labels Character vector of current labels (typically
#'   `object$clusters`).
#' @param mapping Named character vector mapping old labels to new ones.
#' @param keep.original.if.missing If `TRUE` (default), labels not in the
#'   mapping keep their original value. If `FALSE`, they become `NA`.
#'
#' @return Character vector of recoded labels.
#'
#' @examples
#' \dontrun{
#' labels <- c("CD4_T_cells", "CD8_T_cells", "T_reg_cells",
#'             "Neutrophils", "Cancer_CXCL14")
#' mapping <- c(CD4_T_cells = "T_cells",
#'              CD8_T_cells = "T_cells",
#'              T_reg_cells = "T_cells",
#'              Cancer_CXCL14 = "Cancer")
#' CollapseLabels(labels, mapping)
#' # "T_cells" "T_cells" "T_cells" "Neutrophils" "Cancer"
#' }
#'
#' @export
CollapseLabels <- function(labels, mapping, keep.original.if.missing = TRUE) {
  labels <- as.character(labels)
  hits   <- labels %in% names(mapping)
  if (keep.original.if.missing) {
    labels[hits] <- unname(mapping[labels[hits]])
    labels
  } else {
    out <- rep(NA_character_, length(labels))
    out[hits] <- unname(mapping[labels[hits]])
    out
  }
}

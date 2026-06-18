#' Internal utility functions
#'
#' Helpers for layer access (Seurat v4 vs v5), label collapsing, validation,
#' and verbose messaging.
#'
#' @keywords internal
#' @name xeniumclean-utils
NULL

# Get count or data matrix from any Seurat object (v4 or v5 compatible) ------

#' Flexible counts/data accessor
#'
#' Works across Seurat v4 (`slot`) and v5 (`layer`) APIs and handles split layers
#' by joining them on the fly without modifying the original object.
#'
#' @param object A Seurat object.
#' @param assay Assay name. If `NULL`, uses `DefaultAssay`.
#' @param layer Layer name. `"counts"` or `"data"` (or `"scale.data"`).
#' @keywords internal
.get_layer_data <- function(object, assay = NULL, layer = "counts") {
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)

  assay_obj <- object[[assay]]

  # Seurat v5: Assay5 class, may have split layers
  if (methods::is(assay_obj, "Assay5")) {
    available <- SeuratObject::Layers(assay_obj)

    if (layer %in% available) {
      return(SeuratObject::LayerData(object, assay = assay, layer = layer))
    }

    # Try to find split versions like "counts.SAMPLE1"
    split_layers <- grep(paste0("^", layer, "\\."), available, value = TRUE)
    if (length(split_layers) > 0) {
      # Join temporarily to a copy
      tmp <- SeuratObject::JoinLayers(object, assay = assay)
      return(SeuratObject::LayerData(tmp, assay = assay, layer = layer))
    }

    stop("Layer '", layer, "' not found in assay '", assay,
         "'. Available: ", paste(available, collapse = ", "))
  }

  # Seurat v4 Assay
  SeuratObject::GetAssayData(object, slot = layer, assay = assay)
}

# Collapse / recode cluster labels using a named vector ---------------------

#' Apply a named-vector label recoder
#'
#' Returns a character vector the same length as `labels`, with values in
#' `names(mapping)` replaced by the corresponding `mapping` values. Labels not
#' in the mapping keep their original value.
#'
#' @param labels Character vector of current labels.
#' @param mapping Named character vector: `c("old_label" = "new_label")`.
#' @keywords internal
.apply_label_mapping <- function(labels, mapping) {
  labels <- as.character(labels)
  hits <- labels %in% names(mapping)
  labels[hits] <- unname(mapping[labels[hits]])
  labels
}

# Validate a Seurat object has what we need ---------------------------------

#' Validate object setup before cleanup
#'
#' @param object A Seurat object.
#' @param assay Assay name expected to exist.
#' @param group.by Metadata column expected to exist.
#' @keywords internal
.validate_object <- function(object, assay, group.by) {
  if (!methods::is(object, "Seurat")) {
    stop("`object` must be a Seurat object, got ", class(object)[1])
  }
  if (!assay %in% SeuratObject::Assays(object)) {
    stop("Assay '", assay, "' not found. Available: ",
         paste(SeuratObject::Assays(object), collapse = ", "))
  }
  if (!group.by %in% colnames(object@meta.data)) {
    stop("Metadata column '", group.by,
         "' not found. Available: ",
         paste(colnames(object@meta.data), collapse = ", "))
  }
  invisible(TRUE)
}

# Validate the gene sets list -----------------------------------------------

#' Validate a XeniumCleanGeneSets object
#'
#' @param gene.sets A gene sets list returned by [BuildReferenceGeneSets()].
#' @keywords internal
.validate_gene_sets <- function(gene.sets) {
  required <- c("ExpressedGenes", "NotExpressedGenes")
  if (!all(required %in% names(gene.sets))) {
    stop("`gene.sets` must contain at least ExpressedGenes and NotExpressedGenes.")
  }
  if (!is.list(gene.sets$ExpressedGenes) || !is.list(gene.sets$NotExpressedGenes)) {
    stop("ExpressedGenes and NotExpressedGenes must both be lists.")
  }
  invisible(TRUE)
}

# Verbose-aware message wrapper ---------------------------------------------

#' Conditional message
#'
#' @param ... Passed to `message()`.
#' @param verbose Logical. If `FALSE`, suppress.
#' @keywords internal
.v_message <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) message(...)
}

# Skip-pair check -----------------------------------------------------------

#' Are two cell types in the same skip-pair?
#'
#' @param a,b Cell type strings.
#' @param skip.pairs List of length-2 character vectors.
#' @keywords internal
.is_skip_pair <- function(a, b, skip.pairs) {
  if (is.null(skip.pairs) || length(skip.pairs) == 0) return(FALSE)
  any(vapply(skip.pairs, function(p) all(c(a, b) %in% p), logical(1)))
}

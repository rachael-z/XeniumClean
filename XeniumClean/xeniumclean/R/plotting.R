#' Visualise transcripts removed per cell
#'
#' Builds a violin (or jitter) plot of `xeniumclean_removed_transcripts` grouped
#' by cluster, optionally faceted by sample. Useful for spotting clusters where
#' cleanup is acting most aggressively (often a sign of dense neighbour types
#' or potential misclassification).
#'
#' Requires `ggplot2`.
#'
#' @param object A cleaned Seurat object or named list of cleaned objects. Must
#'   contain `xeniumclean_removed_transcripts` (added by [XeniumClean()]).
#' @param group.by Metadata column for the x-axis grouping (cluster column).
#' @param facet.by Optional second metadata column for facet wrapping.
#' @param log.y If `TRUE` (default), log10 y-axis.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' PlotRemovedTranscripts(xenium_list, group.by = "cell_type")
#' }
#'
#' @export
PlotRemovedTranscripts <- function(object,
                                   group.by,
                                   facet.by = NULL,
                                   log.y    = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to use PlotRemovedTranscripts().")
  }

  per_obj <- function(obj, samp = NA) {
    if (!"xeniumclean_removed_transcripts" %in% colnames(obj@meta.data)) {
      stop("Object missing xeniumclean_removed_transcripts; run XeniumClean first.")
    }
    df <- data.frame(
      sample   = samp,
      cluster  = as.character(obj@meta.data[[group.by]]),
      removed  = obj$xeniumclean_removed_transcripts,
      stringsAsFactors = FALSE
    )
    if (!is.null(facet.by) && facet.by %in% colnames(obj@meta.data)) {
      df$facet <- as.character(obj@meta.data[[facet.by]])
    }
    df
  }

  df <- if (methods::is(object, "Seurat")) {
    per_obj(object)
  } else {
    do.call(rbind, lapply(names(object),
                          function(s) per_obj(object[[s]], samp = s)))
  }

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = .data$cluster, y = .data$removed + 1)) +
    ggplot2::geom_violin(scale = "width", fill = "lightgrey") +
    ggplot2::geom_jitter(alpha = 0.05, size = 0.3) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 1)) +
    ggplot2::labs(y = "Transcripts removed per cell (+1)",
                  x = NULL)

  if (isTRUE(log.y)) p <- p + ggplot2::scale_y_log10()

  if (!is.null(facet.by) && "facet" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~ .data$facet, ncol = 3)
  } else if (!methods::is(object, "Seurat")) {
    p <- p + ggplot2::facet_wrap(~ .data$sample, ncol = 3)
  }

  p
}

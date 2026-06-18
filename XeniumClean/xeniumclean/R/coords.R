#' Get spatial coordinates from a Seurat object
#'
#' Flexible coordinate accessor that handles the multiple ways spatial
#' coordinates can be stored in a Seurat object: image FOVs (Xenium standard),
#' a `spatial` dimensional reduction (common for older or non-image-based
#' workflows), or arbitrary metadata columns. Returns a matrix with cells in
#' rows aligned to `colnames(object)`.
#'
#' @param object A Seurat object.
#' @param coords.source One of:
#'   - `"auto"` (default): try image FOVs first, then `spatial` reduction,
#'     then metadata.
#'   - `"image"`: pull from `obj@images[[fov]]@boundaries$centroids`.
#'   - `"reduction"`: pull from `obj@reductions[[reduction.name]]`.
#'   - `"metadata"`: pull from `obj@meta.data[, coord.cols]`.
#'   - `"matrix"`: use the matrix supplied to `coords`.
#' @param fov FOV name when `coords.source = "image"`. Defaults to the first FOV.
#' @param reduction.name Reduction to use when `coords.source = "reduction"`.
#'   Defaults to `"spatial"`.
#' @param coord.cols Length-2 character vector of metadata column names for
#'   `coords.source = "metadata"`. Defaults to `c("x","y")`.
#' @param coords An external numeric matrix with rownames matching
#'   `colnames(object)` and 2 columns (x, y). Used when
#'   `coords.source = "matrix"`.
#'
#' @return A two-column numeric matrix with rownames aligned to
#'   `colnames(object)`.
#'
#' @examples
#' \dontrun{
#' # Standard Xenium object with image FOV
#' coords <- GetSpatialCoords(xenium_obj)
#'
#' # Object where coords live in a 'spatial' reduction
#' coords <- GetSpatialCoords(xenium_obj, coords.source = "reduction",
#'                            reduction.name = "spatial")
#'
#' # Object where coords live in metadata columns named cx, cy
#' coords <- GetSpatialCoords(xenium_obj, coords.source = "metadata",
#'                            coord.cols = c("cx", "cy"))
#' }
#'
#' @export
GetSpatialCoords <- function(object,
                             coords.source  = c("auto","image","reduction","metadata","matrix"),
                             fov            = NULL,
                             reduction.name = "spatial",
                             coord.cols     = c("x", "y"),
                             coords         = NULL) {

  coords.source <- match.arg(coords.source)

  if (coords.source == "matrix") {
    if (is.null(coords)) stop("coords.source = 'matrix' requires `coords`.")
    if (!all(colnames(object) %in% rownames(coords))) {
      stop("Supplied coords matrix does not contain all cells of `object`.")
    }
    return(as.matrix(coords[colnames(object), 1:2, drop = FALSE]))
  }

  # Auto: try each in order
  if (coords.source == "auto") {
    if (length(SeuratObject::Images(object)) > 0) {
      coords.source <- "image"
    } else if (reduction.name %in% names(object@reductions)) {
      coords.source <- "reduction"
    } else if (all(coord.cols %in% colnames(object@meta.data))) {
      coords.source <- "metadata"
    } else {
      stop("Could not auto-detect spatial coordinates. ",
           "Specify `coords.source` explicitly.")
    }
  }

  out <- switch(
    coords.source,
    image = {
      if (is.null(fov)) fov <- SeuratObject::Images(object)[1]
      if (!fov %in% SeuratObject::Images(object)) {
        stop("FOV '", fov, "' not found in object.")
      }
      raw <- SeuratObject::GetTissueCoordinates(object[[fov]][["centroids"]])
      m <- as.matrix(raw[, c("x", "y")])
      rownames(m) <- raw$cell
      m
    },
    reduction = {
      if (!reduction.name %in% names(object@reductions)) {
        stop("Reduction '", reduction.name, "' not found.")
      }
      e <- SeuratObject::Embeddings(object, reduction = reduction.name)
      e[, 1:2, drop = FALSE]
    },
    metadata = {
      if (!all(coord.cols %in% colnames(object@meta.data))) {
        stop("Metadata columns ", paste(coord.cols, collapse = ", "),
             " not found.")
      }
      m <- as.matrix(object@meta.data[, coord.cols, drop = FALSE])
      rownames(m) <- colnames(object)
      m
    }
  )

  # Reorder to object cell order; drop any cells not present in coords
  common <- intersect(colnames(object), rownames(out))
  if (length(common) < ncol(object)) {
    warning("Coordinates missing for ", ncol(object) - length(common),
            " cells. They will be excluded from neighbour computation.")
  }
  out[colnames(object)[colnames(object) %in% common], , drop = FALSE]
}


#' Compute spatial neighbours within a radius
#'
#' Identifies, for each cell, the indices of other cells whose coordinates fall
#' within `radius` (in the same units as the coordinates). Wraps
#' [dbscan::frNN()] and returns a list aligned to the cells of `object`.
#'
#' @param object A Seurat object.
#' @param radius Numeric distance threshold (typically microns for Xenium).
#'   The cleanup defaults assume coordinates in microns; if your coordinates are
#'   in pixels, choose `radius` accordingly.
#' @param coords Optional pre-computed coordinate matrix (skips the auto-detect).
#'   If `NULL`, calls [GetSpatialCoords()] with `...`.
#' @param store.in.misc If `TRUE` (default), save the neighbour list to
#'   `object@misc$xeniumclean_neighbors`.
#' @param ... Passed to [GetSpatialCoords()].
#'
#' @return Either:
#'   - The modified Seurat object with neighbours stored in `@misc`, OR
#'   - A list of integer vectors of length `ncol(object)` (each entry holds the
#'     positional indices of that cell's neighbours, excluding itself), if
#'     `store.in.misc = FALSE`.
#'
#' @examples
#' \dontrun{
#' xenium_obj <- ComputeSpatialNeighbors(xenium_obj, radius = 50)
#' length(xenium_obj@misc$xeniumclean_neighbors)
#' }
#'
#' @export
ComputeSpatialNeighbors <- function(object,
                                    radius        = 50,
                                    coords        = NULL,
                                    store.in.misc = TRUE,
                                    ...) {

  if (is.null(coords)) {
    coords <- GetSpatialCoords(object, ...)
  }
  if (nrow(coords) != ncol(object)) {
    stop("Coordinate matrix has ", nrow(coords),
         " rows but object has ", ncol(object), " cells.")
  }

  fr <- dbscan::frNN(coords, eps = radius, sort = FALSE)
  neighbors_list <- lapply(seq_along(fr$id),
                           function(i) setdiff(fr$id[[i]], i))
  neighbors_list <- lapply(neighbors_list, as.integer)

  if (isTRUE(store.in.misc)) {
    object@misc$xeniumclean_neighbors <- neighbors_list
    object@misc$xeniumclean_neighbors_radius <- radius
    return(object)
  }
  neighbors_list
}

#' Spatial neighbour-aware transcript cleanup
#'
#' Removes biologically implausible transcripts from a spatial transcriptomics
#' object. For each cell of type A, identifies its spatial neighbours of other
#' types B, and erases transcripts of genes that are expressed in B but not
#' expressed in A according to the reference gene sets. Genuine A-expressed
#' transcripts are never removed.
#'
#' Works on a single Seurat object (one tissue section) or a list of objects
#' (multiple sections). Cleanup is always performed per-section because spatial
#' neighbours are within-section.
#'
#' @section Before running, confirm:
#' \itemize{
#'   \item \strong{Cell types are pre-defined} in your spatial object via the
#'     `group.by` column. XeniumClean does not cluster or annotate.
#'   \item \strong{Cell type names match} between your spatial data and
#'     `gene.sets`. Use [CollapseLabels()] to align.
#'   \item \strong{Cell types are BROAD} (T_cells, B_cells, Cancer), not fine
#'     subsets (CD4_T_cells, Memory_B_cells, Cancer_Basal). Fine subclustering
#'     should happen AFTER cleanup on the cleaned assay.
#'   \item \strong{Multi-sample objects are split into a list}. Use
#'     `SplitObject(obj, split.by = "sample")` rather than `subset()` (the
#'     latter can be slow and produce inconsistent results on v5 objects with
#'     split layers).
#'   \item \strong{Your reference matches your tissue/disease}. A reference
#'     mismatch is the most common cause of biology-erasing over-cleanup.
#' }
#'
#' @section Subclustering benefit:
#' A key reason to use XeniumClean is that subclustering on the cleaned data
#' reveals real biology. In raw spatial data, fine subclusters often reflect
#' neighbourhood contamination (T cells near macrophages cluster together
#' because they share macrophage transcripts, not because they're a true
#' subtype). After cleanup, subclustering recovers genuine biological subsets
#' (CD4 vs CD8, naive vs effector, etc.).
#'
#' @param object A Seurat object **or a named list of Seurat objects**, one per
#'   tissue section.
#' @param gene.sets A `XeniumCleanGeneSets` object from
#'   [BuildReferenceGeneSets()], or a manually constructed list with
#'   `ExpressedGenes` and `NotExpressedGenes` elements.
#' @param group.by Metadata column on `object` holding cell type labels. These
#'   should map to the names in `gene.sets$ExpressedGenes`. Cell types not in
#'   the gene sets pass through untouched (safe for unmappable populations like
#'   Neutrophils when not in the reference).
#' @param assay Source assay holding the raw spatial counts. Default `"Xenium"`.
#' @param layer Layer/slot holding raw counts. Default `"counts"`.
#' @param new.assay.name Name for the new cleaned assay. Default `"XeniumClean"`.
#' @param radius Spatial radius for neighbour definition, in coordinate units
#'   (microns for Xenium standard). Default `50`.
#' @param coords.source,fov,reduction.name,coord.cols Passed to
#'   [GetSpatialCoords()] when coordinates need to be located. See that function
#'   for details.
#' @param skip.pairs Optional list of length-2 character vectors specifying
#'   cell type pairs to skip cleanup between (useful for closely-related
#'   lineages where one's "expressed" markers really are shared, e.g.
#'   `list(c("T_cells","NK_cells"), c("B_cells","Plasma_cells"))`).
#' @param block.size Cells per parallel block. Default `500`.
#' @param parallel Use `future` for parallel execution. Default `FALSE`.
#' @param workers Number of parallel workers if `parallel = TRUE`. Default `4`.
#' @param verbose Logical; print progress. Default `TRUE`.
#'
#' @return The input object with a new assay (`new.assay.name`) added and a
#'   metadata column `xeniumclean_removed_transcripts` recording how many
#'   transcripts were erased per cell. For a list input, returns a list with
#'   each element cleaned.
#'
#' @section Algorithm summary:
#' For each cell `i` of type A in section `s`:
#' 1. Find spatial neighbours within `radius` of cell `i`.
#' 2. Determine the unique cell types of those neighbours, excluding A itself
#'    and any skip-paired types.
#' 3. For each neighbour cell type B, compute the "forbidden set":
#'    `intersect(ExpressedGenes[[B]], NotExpressedGenes[[A]])`. These are
#'    genes that could plausibly have bled over from a B neighbour but cannot
#'    originate from A itself.
#' 4. Zero out those genes in cell `i`'s count vector.
#'
#' Cells whose label is not in the gene sets are skipped (their counts are
#' preserved). The cleanup is conservative: a gene is only removed if there is
#' a specific nearby cell type that genuinely expresses it AND the focal cell's
#' own type genuinely doesn't.
#'
#' @section Cell name and coordinate alignment:
#' This function does not modify cell ordering. Coordinates are reordered
#' internally to match `colnames(object)`. If any cells have no coordinates,
#' they are skipped (their counts preserved).
#'
#' @examples
#' \dontrun{
#' # Single section
#' xenium_obj <- XeniumClean(
#'   object    = xenium_obj,
#'   gene.sets = gene_sets,
#'   group.by  = "cell_type",
#'   radius    = 50
#' )
#'
#' # Multi-section
#' xenium_list <- XeniumClean(
#'   object    = xenium_list,
#'   gene.sets = gene_sets,
#'   group.by  = "cell_type",
#'   radius    = 50,
#'   parallel  = TRUE,
#'   workers   = 6
#' )
#'
#' # With skip pairs for lineage-related cell types
#' xenium_obj <- XeniumClean(
#'   object    = xenium_obj,
#'   gene.sets = gene_sets,
#'   group.by  = "cell_type",
#'   skip.pairs = list(
#'     c("T_cells",   "NK_cells"),
#'     c("B_cells",   "Plasma_cells"),
#'     c("Monocytes", "Macrophage")
#'   )
#' )
#' }
#'
#' @export
XeniumClean <- function(object,
                        gene.sets,
                        group.by,
                        assay           = "Xenium",
                        layer           = "counts",
                        new.assay.name  = "XeniumClean",
                        radius          = 50,
                        coords.source   = "auto",
                        fov             = NULL,
                        reduction.name  = "spatial",
                        coord.cols      = c("x", "y"),
                        skip.pairs      = NULL,
                        block.size      = 500,
                        parallel        = FALSE,
                        workers         = 4,
                        verbose         = TRUE) {

  .validate_gene_sets(gene.sets)

  # Dispatch on input type --------------------------------------------------

  if (methods::is(object, "Seurat")) {
    return(.clean_single(
      object         = object,
      gene.sets      = gene.sets,
      group.by       = group.by,
      assay          = assay,
      layer          = layer,
      new.assay.name = new.assay.name,
      radius         = radius,
      coords.source  = coords.source,
      fov            = fov,
      reduction.name = reduction.name,
      coord.cols     = coord.cols,
      skip.pairs     = skip.pairs,
      block.size     = block.size,
      parallel       = parallel,
      workers        = workers,
      verbose        = verbose
    ))
  }

  if (is.list(object)) {
    if (is.null(names(object))) {
      stop("Input list must be named (one name per sample).")
    }
    out <- object
    for (samp in names(object)) {
      .v_message("\n=== Cleaning ", samp, " ===", verbose = verbose)
      out[[samp]] <- .clean_single(
        object         = object[[samp]],
        gene.sets      = gene.sets,
        group.by       = group.by,
        assay          = assay,
        layer          = layer,
        new.assay.name = new.assay.name,
        radius         = radius,
        coords.source  = coords.source,
        fov            = fov,
        reduction.name = reduction.name,
        coord.cols     = coord.cols,
        skip.pairs     = skip.pairs,
        block.size     = block.size,
        parallel       = parallel,
        workers        = workers,
        verbose        = verbose
      )
    }
    return(out)
  }

  stop("`object` must be a Seurat object or a named list of Seurat objects.")
}


# Internal: cleanup for one Seurat object -----------------------------------

.clean_single <- function(object,
                          gene.sets,
                          group.by,
                          assay,
                          layer,
                          new.assay.name,
                          radius,
                          coords.source,
                          fov,
                          reduction.name,
                          coord.cols,
                          skip.pairs,
                          block.size,
                          parallel,
                          workers,
                          verbose) {

  .validate_object(object, assay, group.by)

  # 1. Coordinates and neighbours
  coords <- GetSpatialCoords(object,
                             coords.source  = coords.source,
                             fov            = fov,
                             reduction.name = reduction.name,
                             coord.cols     = coord.cols)
  # Restrict object to cells with coords (rare edge case)
  if (nrow(coords) < ncol(object)) {
    keep_cells <- intersect(colnames(object), rownames(coords))
    object     <- subset(object, cells = keep_cells)
  }
  coords <- coords[colnames(object), , drop = FALSE]

  .v_message("  computing neighbours for ", nrow(coords),
             " cells (radius = ", radius, ")...", verbose = verbose)
  fr <- dbscan::frNN(coords, eps = radius, sort = FALSE)
  neighbors_list <- lapply(seq_along(fr$id),
                           function(i) as.integer(setdiff(fr$id[[i]], i)))

  # 2. Pull source counts (handles v4/v5)
  expr_mat <- .get_layer_data(object, assay = assay, layer = layer)
  celltypes_vec <- as.character(object@meta.data[[group.by]])

  ExpressedGenes    <- gene.sets$ExpressedGenes
  NotExpressedGenes <- gene.sets$NotExpressedGenes
  gene_universe     <- rownames(expr_mat)

  # 3. Cleanup loop, optionally parallel
  n_cells <- ncol(expr_mat)
  blocks  <- split(seq_len(n_cells),
                   ceiling(seq_len(n_cells)/block.size))
  .v_message("  cleaning ", n_cells, " cells in ", length(blocks),
             " blocks of ~", block.size, "...", verbose = verbose)

  worker <- function(block_cells) {
    block_mat <- expr_mat[, block_cells, drop = FALSE]

    for (i in seq_along(block_cells)) {
      col_idx   <- block_cells[i]
      cell_type <- celltypes_vec[col_idx]
      neigh_idx <- neighbors_list[[col_idx]]
      if (length(neigh_idx) == 0) next
      if (!cell_type %in% names(NotExpressedGenes)) next

      neigh_types <- setdiff(unique(celltypes_vec[neigh_idx]), cell_type)
      if (length(neigh_types) == 0) next

      genes_to_remove <- character(0)
      for (nt in neigh_types) {
        if (!nt %in% names(ExpressedGenes)) next
        if (.is_skip_pair(cell_type, nt, skip.pairs)) next
        forbidden <- intersect(ExpressedGenes[[nt]],
                               NotExpressedGenes[[cell_type]])
        genes_to_remove <- union(genes_to_remove, forbidden)
      }

      if (length(genes_to_remove) == 0) next
      genes_to_remove <- intersect(genes_to_remove, gene_universe)
      if (length(genes_to_remove) > 0) {
        block_mat[genes_to_remove, i] <- 0
      }
    }
    block_mat
  }

  blocks_out <- if (isTRUE(parallel) &&
                   requireNamespace("future.apply", quietly = TRUE)) {
    if (requireNamespace("future", quietly = TRUE)) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = workers)
    }
    future.apply::future_lapply(blocks, worker, future.seed = TRUE)
  } else {
    lapply(blocks, worker)
  }

  expr_mat_clean <- do.call(cbind, blocks_out)
  dimnames(expr_mat_clean) <- dimnames(expr_mat)

  # 4. Attach as a new assay and bookkeeping
  object[[new.assay.name]] <- SeuratObject::CreateAssayObject(
                                counts = expr_mat_clean)

  object$xeniumclean_removed_transcripts <- as.numeric(
    Matrix::colSums(expr_mat - expr_mat_clean))

  pct <- 100 * sum(object$xeniumclean_removed_transcripts) / sum(expr_mat)
  .v_message(sprintf("  removed %s transcripts (%.2f%% of total)",
                     format(sum(object$xeniumclean_removed_transcripts),
                            big.mark = ","),
                     pct),
             verbose = verbose)

  object
}


#' Summarise cleanup results
#'
#' Builds a data frame describing how many transcripts were removed per sample
#' and per cluster, plus rates.
#'
#' @param object A cleaned Seurat object or list of cleaned objects.
#' @param group.by Metadata column for the cluster split. Default
#'   `"xclean_label"` (a common convention) but pass whatever you used for
#'   [XeniumClean()].
#' @param source.assay Name of the original (uncleaned) assay. Default
#'   `"Xenium"`.
#' @param cleaned.assay Name of the cleaned assay. Default `"XeniumClean"`.
#'
#' @return A data frame with one row per sample (or per `sample_x_cluster`
#'   combination) summarising removed-transcript counts and percentages.
#'
#' @examples
#' \dontrun{
#' summary_df <- CleanupSummary(xenium_list,
#'                              group.by = "cell_type",
#'                              source.assay = "Xenium",
#'                              cleaned.assay = "XeniumClean")
#' summary_df
#' }
#'
#' @export
CleanupSummary <- function(object,
                           group.by      = NULL,
                           source.assay  = "Xenium",
                           cleaned.assay = "XeniumClean") {

  per_sample <- function(obj, samp = NA_character_) {
    orig  <- .get_layer_data(obj, assay = source.assay,  layer = "counts")
    clean <- .get_layer_data(obj, assay = cleaned.assay, layer = "counts")
    df <- data.frame(
      sample          = samp,
      total_cells     = ncol(obj),
      total_reads     = sum(orig),
      reads_removed   = sum(orig - clean),
      pct_removed     = round(100 * sum(orig - clean) / sum(orig), 2),
      median_per_cell = stats::median(Matrix::colSums(orig - clean)),
      max_per_cell    = max(Matrix::colSums(orig - clean))
    )
    if (!is.null(group.by) && group.by %in% colnames(obj@meta.data)) {
      per_clu <- tapply(
        Matrix::colSums(orig - clean),
        as.character(obj@meta.data[[group.by]]),
        function(x) c(median = stats::median(x), max = max(x), n = length(x))
      )
      df$per_cluster <- list(do.call(rbind, per_clu))
    }
    df
  }

  if (methods::is(object, "Seurat")) return(per_sample(object))

  do.call(rbind, lapply(names(object),
                        function(s) per_sample(object[[s]], samp = s)))
}

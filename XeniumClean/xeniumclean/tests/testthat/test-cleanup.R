test_that("CollapseLabels recodes correctly", {
  lab <- c("CD4_T_cells","CD8_T_cells","B_cells","Cancer_LumA")
  map <- c(CD4_T_cells = "T_cells",
           CD8_T_cells = "T_cells",
           Cancer_LumA = "Cancer")
  out <- CollapseLabels(lab, map)
  expect_equal(out, c("T_cells","T_cells","B_cells","Cancer"))
})

test_that("CollapseLabels keeps missing labels when requested", {
  lab <- c("X","Y","Z")
  map <- c(X = "A")
  expect_equal(CollapseLabels(lab, map, keep.original.if.missing = TRUE),
               c("A","Y","Z"))
  expect_equal(CollapseLabels(lab, map, keep.original.if.missing = FALSE),
               c("A", NA, NA))
})

test_that("BuildReferenceGeneSets fails on bad thresholds", {
  skip_if_not_installed("SeuratObject")
  fake_obj <- structure(list(),
                        class = "Seurat",
                        meta.data = data.frame(group = factor(c("a","b"))))
  expect_error(
    BuildReferenceGeneSets(fake_obj, group.by = "group",
                            expressed.threshold     = 0.05,
                            not.expressed.threshold = 0.10),
    "must be greater"
  )
})

test_that("XeniumClean errors gracefully on bad gene.sets", {
  expect_error(
    XeniumClean(object = "not_a_seurat",
                gene.sets = list(),
                group.by = "x"),
    "ExpressedGenes"
  )
})

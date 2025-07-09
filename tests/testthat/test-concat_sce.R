
test_concat_names <- function(obj.list, merged.obj) {
  expect_equal(ncol(merged.obj), sum(sapply(obj.list, ncol)))
  expect_setequal(
    rownames(merged.obj),
    unique(unlist(lapply(obj.list, rownames)))
  )
  expect_equal(anyDuplicated(rownames(merged.obj)), 0)
  expect_equal(anyDuplicated(colnames(merged.obj)), 0)

  expect_setequal(
    colnames(rowData(merged.obj)),
    unique(unlist(lapply(obj.list, function(x) colnames(rowData(x)))))
  )
  expect_setequal(
    colnames(colData(merged.obj)),
    unique(unlist(lapply(obj.list, function(x) colnames(colData(x)))))
  )
  invisible(NULL)
}

test_concat_add_idx <- function(obj.list, merged.obj, add.idx) {
  expect_setequal(
    colnames(merged.obj),
    unique(unlist(lapply(add.idx, function(i) {
      paste0(i, "_", colnames(obj.list[[i]]))
    })))
  )
  invisible(NULL)
}

test_concat_assays <- function(obj.list, merged.obj) {
  expect_setequal(
    assayNames(merged.obj),
    Reduce(intersect, lapply(obj.list, assayNames))
  )
  for (i in assayNames(merged.obj)) {
    m0 <- assay(merged.obj, i = i)
    m.list <- lapply(obj.list, assay, i = i)
    m1 <- concatByCols(m.list[[1]], m.list[-1])
    expect_equal(unname(m1[rownames(m0), ]), unname(m0))
  }
  invisible(NULL)
}

test_concat_reds <- function(obj.list, merged.obj) {
  expect_setequal(
    reducedDimNames(merged.obj),
    unlist(Reduce(intersect, lapply(obj.list, reducedDimNames)))
  )
}

test_concat_FE <- function(obj.list) {
  merged.obj <- concatByCols(obj.list[[1]], obj.list[-1], verbose = FALSE)
  test_concat_names(obj.list, merged.obj)
  test_concat_assays(obj.list, merged.obj)

  add.idx <- names(obj.list)
  merged.obj <- concatByCols(
    obj.list[[1]], obj.list[-1],
    add.idx = add.idx,
    verbose = FALSE
  )
  test_concat_add_idx(obj.list, merged.obj, add.idx)
  test_concat_assays(obj.list, merged.obj)
  merged.obj
}

test_concat_SCFE <- function(obj.list) {
  merged.obj <- test_concat_FE(obj.list)
  test_concat_reds(obj.list, merged.obj)
  merged.obj
}

test_that("concatByCols for FlexExperiment", {
  sces <- sapply(.all_examples, exampleFE)
  m <- test_concat_FE(sces)

  # Set duplicated cell names
  colnames(sces[[2]])[1:10] <- colnames(sces[[3]])[1:10]
  expect_warning(m <- test_concat_FE(sces))

  # Set missing colData
  colData(sces[[1]])$seurat_clusters <- NULL
  expect_warning(m <- test_concat_FE(sces))
  expect_true(anyNA(colData(m)$seurat_clusters))

  # Set missing rowData
  rowData(sces[[1]])$Type <- NULL
  expect_warning(m <- test_concat_FE(sces))
  expect_true(anyNA(rowData(m)$Type))
})

test_that("concatByCols for SummarizedExperiment", {
  sces <- sapply(.all_examples, exampleSE)
  m <- test_concat_FE(sces)

  # Set duplicated cell names
  colnames(sces[[2]])[1:10] <- colnames(sces[[3]])[1:10]
  expect_warning(m <- test_concat_FE(sces))

  # Set missing colData
  colData(sces[[1]])$seurat_clusters <- NULL
  expect_warning(m <- test_concat_FE(sces))
  expect_true(anyNA(colData(m)$seurat_clusters))

  # Set missing rowData
  rowData(sces[[1]])$Type <- NULL
  expect_warning(m <- test_concat_FE(sces))
  expect_true(anyNA(rowData(m)$Type))
})

test_that("concatByCols for SingleCellExperiment", {
  sces <- sapply(.all_examples, exampleSCE)
  for (i in seq_along(sces)) {
    altExp(sces[[i]], "raw") <- sces[[i]][1:10, ]
  }
  m <- test_concat_SCFE(sces)

  expect_equal(
    altExpNames(m),
    unlist(Reduce(intersect, lapply(sces, altExpNames)))
  )
  expect_equal(colnames(altExp(m)), colnames(m))

  # Set duplicated cell names
  colnames(sces[[2]])[1:10] <- colnames(sces[[3]])[1:10]
  expect_warning(m <- test_concat_SCFE(sces))

  # Set missing colData
  colData(sces[[1]])$seurat_clusters <- NULL
  expect_warning(m <- test_concat_SCFE(sces))
  expect_true(anyNA(colData(m)$seurat_clusters))

  # Set missing rowData
  rowData(sces[[1]])$Type <- NULL
  expect_warning(m <- test_concat_SCFE(sces))
  expect_true(anyNA(rowData(m)$Type))

  # Set missing reducedDims
  reducedDim(sces[[1]]) <- NULL
  expect_warning(m <- test_concat_SCFE(sces))
  expect_equal(
    altExpNames(m),
    unlist(Reduce(intersect, lapply(sces, altExpNames)))
  )
  expect_equal(colnames(altExp(m)), colnames(m))
})

test_that("concatByCols for SingleCellFlexExperiment", {
  sces <- sapply(.all_examples, exampleSCFE)
  m <- test_concat_SCFE(sces)

  expect_setequal(
    reducedDimNames(m),
    unlist(Reduce(intersect, lapply(sces, reducedDimNames)))
  )

  # Set duplicated cell names
  colnames(sces[[2]])[1:10] <- colnames(sces[[3]])[1:10]
  expect_warning(m <- test_concat_SCFE(sces))

  # Set missing colData
  colData(sces[[1]])$seurat_clusters <- NULL
  expect_warning(m <- test_concat_SCFE(sces))
  expect_true(anyNA(colData(m)$seurat_clusters))

  # Set missing rowData
  rowData(sces[[1]])$Type <- NULL
  expect_warning(m <- test_concat_SCFE(sces))
  expect_true(anyNA(rowData(m)$Type))
})



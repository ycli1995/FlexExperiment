
add_rowPair <- function(sce) {
  cor1 <- Matrix::rsparsematrix(nrow(sce), nrow(sce), density = 0.3)
  cor2 <- Matrix::rsparsematrix(nrow(sce), nrow(sce), density = 0.3)
  colnames(cor1) <- colnames(cor2) <- rownames(sce)
  rownames(cor1) <- rownames(cor2) <- rownames(sce)
  rowPair(sce, "cor1") <- cor1
  rowPair(sce, "cor2") <- cor2
  sce
}

get_int_test <- function(sce, getfun, xfun, vfun, class = NULL, ...) {
  out <- getfun(sce, ...)
  if (!is.null(class)) {
    expect_true(inherits(out, class))
  }
  expect_equal_no_attr(vfun(out), xfun(sce))
  return(invisible(NULL))
}

get_colPair_test <- function(sce, ...) {
  class <- "SelfHits"
  get_int_test(sce, colPair, ncol, nnode, class, ...)

  class <- "CsparseMatrix"
  get_int_test(sce, colPair, colnames, colnames, class, asSparse = TRUE)
  get_int_test(sce, colPair, colnames, rownames, class, asSparse = TRUE)
}

get_rowPair_test <- function(sce, ...) {
  class <- "SelfHits"
  get_int_test(sce, rowPair, nrow, nnode, class, ...)

  class <- "CsparseMatrix"
  get_int_test(sce, rowPair, rownames, colnames, class, asSparse = TRUE)
  get_int_test(sce, rowPair, rownames, rownames, class, asSparse = TRUE)
}

test_subset_SCE <- function(obj, i = NULL, j = NULL) {
  reducs <- reducedDims(obj)
  colps <- colPairs(obj, asSparse = TRUE)
  rowps <- rowPairs(obj, asSparse = TRUE)

  rr <- rowRanges(obj)
  rd <- rowData(obj)
  cd <- colData(obj)

  if (is.null(i)) {
    if (!is.null(j)) {
      obj2 <- obj[, j]
    }
  } else {
    if (is.null(j)) {
      obj2 <- obj[i, ]
    } else {
      obj2 <- obj[i, j]
    }
  }

  rr2 <- rowRanges(obj2)
  rd2 <- rowData(obj2)
  cd2 <- colData(obj2)

  reducs2 <- reducedDims(obj2)
  colps2 <- colPairs(obj2, asSparse = TRUE)
  rowps2 <- rowPairs(obj2, asSparse = TRUE)

  if (!is.null(i)) {
    expect_identical(rr[i], rr2)
    expect_identical(rd[i, , drop = FALSE], rd2)
    for (ii in seq_along(rowps)) {
      expect_identical(rowps2[[ii]], rowps[[ii]][i, i, drop = FALSE])
    }
  }
  if (!is.null(j)) {
    expect_identical(cd[j, , drop = FALSE], cd2)
    for (ii in seq_along(colps)) {
      expect_identical(colps2[[ii]], colps[[ii]][j, j, drop = FALSE])
    }
    for (ii in seq_along(reducs)) {
      expect_identical(reducs2[[ii]], reducs[[ii]][j, , drop = FALSE])
    }
  }
  return(obj2)
}

test_that("reducedDims getters", {
  sce <- exampleSCFE()

  class <- "LinearEmbeddingMatrix"
  pca <- get_int_test(sce, reducedDim, colnames, rownames, class)
  pca <- get_int_test(sce, reducedDim, colnames, rownames, class, type = "PCA")
  pca <- get_int_test(sce, reducedDim, colnames, rownames, class, type = 1)
  pca <- get_int_test(sce, reducedDim, colnames, rownames, class, type = 1L)

  class <- "matrix"
  umap <- get_int_test(sce, reducedDim, colnames, rownames, class, "UMAP")
  umap <- get_int_test(sce, reducedDim, colnames, rownames, class, 3)
  umap <- get_int_test(sce, reducedDim, colnames, rownames, class, 3L)

  all.reds <- reducedDims(sce)
  expect_s4_class(all.reds, "SimpleList")
  expect_equal_no_attr(names(all.reds), reducedDimNames(sce))
  expect_equal_no_attr(
    as.list(all.reds),
    lapply(reducedDimNames(sce), reducedDim, x = sce)
  )

  expect_equal_no_attr(reducedDim(sce), reducedDims(sce)[[1]])
  expect_equal_no_attr(reducedDim(sce, 2), reducedDims(sce)[["TSNE"]])
  expect_equal_no_attr(reducedDim(sce, "UMAP"), reducedDims(sce)[[3]])
})

test_that("reducedDims<- set NULL", {
  sce <- exampleSCFE()

  reducedDims(sce) <- NULL
  expect_in("reducedDims", colnames(int_colData(sce)))
  expect_length(reducedDims(sce), 0)
  expect_length(reducedDimNames(sce), 0)
})

test_that("reducedDims<- set new list", {
  sce <- exampleSCFE()

  old.reds <- reducedDims(sce)
  reducedDims(sce) <- old.reds
  expect_equal_no_attr(reducedDims(sce), old.reds)

  ## new reducedDims without names
  new.reds <- unname(old.reds)
  expect_warning(reducedDims(sce) <- new.reds, "unnamed")
  expect_equal(reducedDimNames(sce), paste0("unnamed", seq_along(new.reds)))

  new.reds <- list(PCA = old.reds[["PCA"]], old.reds[[2]], UMAP = old.reds[[3]])
  expect_warning(reducedDims(sce) <- new.reds, "unnamed")
  expect_equal(reducedDimNames(sce)[2], "unnamed2")
  expect_equal(reducedDimNames(sce)[c(1, 3)], names(old.reds)[c(1, 3)])

  ## new reducedDims without rownames
  new.reds <- lapply(old.reds, unname)
  expect_error(reducedDims(sce) <- new.reds)

  new.reds <- old.reds
  new.reds[[2]] <- unname(new.reds[[2]])
  expect_error(reducedDims(sce) <- new.reds)

  ## new reducedDims when 'withDimnames = FALSE'
  new.reds <- lapply(old.reds, unname)
  reducedDims(sce, withDimnames = FALSE) <- new.reds
  expect_equal(lapply(reducedDims(sce), rownames), lapply(old.reds, rownames))

  new.reds <- old.reds
  new.reds[[2]] <- unname(new.reds[[2]])
  reducedDims(sce, withDimnames = FALSE) <- new.reds
  expect_equal(lapply(reducedDims(sce), rownames), lapply(old.reds, rownames))

  ## new reducedDims with wrong nrow
  new.reds <- old.reds
  new.reds[[2]] <- head(new.reds[[2]])
  expect_error(reducedDims(sce, withDimnames = FALSE) <- new.reds)

  ## reducedDims(sce)[[2]] <-
  expect_warning(reducedDims(sce)[[4]] <- old.reds[[1]])
  expect_equal_no_attr(reducedDims(sce)[[4]], old.reds[[1]])

  reducedDims(sce)[["PCA2"]] <- old.reds[[1]]
  expect_equal_no_attr(reducedDims(sce)[["PCA2"]], old.reds[[1]])
})

test_that("reducedDim<- type = 'missing'", {
  sce <- exampleSCFE()

  old.reds <- reducedDims(sce)

  ## Modify an existing reducedDim
  reducedDim(sce) <- old.reds[[2]]
  expect_equal(reducedDim(sce), old.reds[[2]])
  expect_equal(reducedDim(sce, 1), old.reds[[2]])
  expect_equal(reducedDim(sce, reducedDimNames(sce)[[1]]), old.reds[[2]])

  ## `withDimnames = TRUE` works
  new.red <- unname(old.reds[[2]])
  expect_error(reducedDim(sce) <- new.red)

  ## `withDimnames = FALSE` works
  reducedDim(sce, withDimnames = FALSE) <- new.red
  expect_equal(rownames(reducedDim(sce)), colnames(sce))
  expect_equal(reducedDim(sce, withDimnames = FALSE), new.red)

  ## remove all reducedDims first
  reducedDims(sce) <- NULL
  reducedDim(sce) <- old.reds[[2]]
  expect_equal(reducedDim(sce), old.reds[[2]])
  expect_equal(reducedDim(sce, 1), old.reds[[2]])
  expect_equal(reducedDim(sce, "unnamed1"), old.reds[[2]])
  expect_equal(reducedDimNames(sce), "unnamed1")
})

test_that("reducedDim<- type = 'numeric'", {
  sce <- exampleSCFE()

  old.reds <- reducedDims(sce)

  ## Modify an existing reducedDim
  reducedDim(sce, 3) <- old.reds[[2]]
  expect_equal(reducedDim(sce, 3), old.reds[[2]])
  expect_equal(reducedDim(sce, reducedDimNames(sce)[[3]]), old.reds[[2]])
  expect_equal(reducedDimNames(sce), names(old.reds))

  ## withDimnames = TRUE
  new.red <- unname(old.reds[[2]])
  expect_error(reducedDim(sce, 3) <- new.red)

  ## withDimnames = FALSE
  reducedDim(sce, 3, withDimnames = FALSE) <- new.red
  expect_equal(rownames(reducedDim(sce, 3)), colnames(sce))
  expect_equal(reducedDim(sce, 3, withDimnames = FALSE), new.red)

  ## Cannot add a new reducedDim with numeric index
  expect_error(reducedDim(sce, 4) <- old.reds[[2]])

  ## Remove all reducedDims first. Still cannot add a new reducedDim
  reducedDims(sce) <- NULL
  expect_error(reducedDim(sce, 1) <- old.reds[[2]])
})

test_that("reducedDim<- type = 'character'", {
  sce <- exampleSCFE()

  old.reds <- reducedDims(sce)

  ## Modify an existing reducedDim
  reducedDim(sce, reducedDimNames(sce)[[3]]) <- old.reds[[2]]
  expect_equal(reducedDim(sce, 3), old.reds[[2]])
  expect_equal(reducedDim(sce, reducedDimNames(sce)[[3]]), old.reds[[2]])
  expect_equal(reducedDimNames(sce), names(old.reds))

  ## withDimnames = TRUE
  new.red <- unname(old.reds[[2]])
  expect_error(reducedDim(sce, reducedDimNames(sce)[[3]]) <- new.red)

  ## withDimnames = FALSE
  reducedDim(sce, reducedDimNames(sce)[[3]], withDimnames = FALSE) <- new.red
  expect_equal(rownames(reducedDim(sce, 3)), colnames(sce))
  expect_equal(reducedDim(sce, 3, withDimnames = FALSE), new.red)

  ## Can add a new reducedDim with character index
  reducedDim(sce, "new") <- old.reds[[2]]
  expect_equal(reducedDim(sce, 4), old.reds[[2]])
  expect_equal(reducedDim(sce, reducedDimNames(sce)[[4]]), old.reds[[2]])
  expect_equal(reducedDimNames(sce), c(names(old.reds), "new"))

  ## Remove all reducedDims, then add a new reducedDim
  reducedDims(sce) <- NULL
  reducedDim(sce, "new") <- old.reds[[2]]
  expect_equal(reducedDim(sce), old.reds[[2]])
  expect_equal(reducedDim(sce, 1), old.reds[[2]])
  expect_equal(reducedDim(sce, reducedDimNames(sce)[[1]]), old.reds[[2]])
  expect_equal(reducedDimNames(sce), "new")
})

test_that("reducedDimNames<-", {
  sce <- exampleSCFE()

  old.reds <- reducedDims(sce)

  ## Modify an existing reducedDimNames
  reducedDimNames(sce)[2] <- "tsne"
  expect_equal(reducedDimNames(sce)[2], "tsne")
  expect_equal(reducedDimNames(sce)[c(1, 3)], names(old.reds)[c(1, 3)])

  ## Remove all reducedDims
  expect_warning(reducedDimNames(sce) <- character())
  expect_equal(reducedDimNames(sce), paste0("unnamed", seq(reducedDims(sce))))
  expect_error(reducedDimNames(sce) <- NULL)
})

test_that("colPairs getters", {
  sce <- exampleSCFE()

  get_colPair_test(sce)
  get_colPair_test(sce, 1)
  get_colPair_test(sce, 1L)
  get_colPair_test(sce, "KNN")
  get_colPair_test(sce, 2)
  get_colPair_test(sce, 2L)
  get_colPair_test(sce, "SNN")

  all.cp <- colPairs(sce)
  expect_s4_class(all.cp, "SimpleList")
  expect_equal_no_attr(names(all.cp), colPairNames(sce))
  expect_equal_no_attr(
    as.list(all.cp),
    lapply(colPairNames(sce), colPair, x = sce)
  )

  ## `asSparse = TRUE` works
  all.cp <- colPairs(sce, asSparse = TRUE)
  expect_s4_class(all.cp, "SimpleList")
  expect_equal_no_attr(names(all.cp), colPairNames(sce))
  expect_equal_no_attr(
    as.list(all.cp),
    lapply(colPairNames(sce), colPair, x = sce, asSparse = TRUE)
  )

  expect_equal_no_attr(colPair(sce), colPairs(sce)[[1]])
  expect_equal_no_attr(colPair(sce, 2), colPairs(sce)[["SNN"]])
  expect_equal_no_attr(colPair(sce, "SNN"), colPairs(sce)[[2]])

  expect_equal_no_attr(
    colPair(sce, asSparse = TRUE),
    colPairs(sce, asSparse = TRUE)[[1]]
  )
  expect_equal_no_attr(
    colPair(sce, 2, asSparse = TRUE),
    colPairs(sce, asSparse = TRUE)[["SNN"]]
  )
  expect_equal_no_attr(
    colPair(sce, "SNN", asSparse = TRUE),
    colPairs(sce, asSparse = TRUE)[[2]]
  )
})

test_that("colPairs<- set NULL", {
  sce <- exampleSCFE()

  colPairs(sce) <- NULL
  expect_in("colPairs", colnames(int_colData(sce)))
  expect_length(colPairs(sce), 0)
  expect_length(colPairNames(sce), 0)
})

test_that("colPairs<- set new list", {
  sce <- exampleSCFE()

  old.cp <- colPairs(sce)
  colPairs(sce) <- old.cp
  expect_equal_no_attr(colPairs(sce), old.cp)

  ## new colPairs without names
  new.cp <- unname(old.cp)
  expect_warning(colPairs(sce) <- new.cp, "unnamed")
  expect_equal(colPairNames(sce), paste0("unnamed", seq(new.cp)))

  new.cp <- list(KNN = old.cp[["KNN"]], old.cp[[2]])
  expect_warning(colPairs(sce) <- new.cp, "unnamed")
  expect_equal(colPairNames(sce)[2], "unnamed2")
  expect_equal(colPairNames(sce)[1], names(old.cp)[1])

  ## new colPairs with wrong nnode
  new.cp <- old.cp
  new.cp[[1]] <- PairedHits(new.cp[[1]])[1:10]@hits
  expect_error(colPairs(sce) <- new.cp)

  ## new colPairs as sparseMatrix
  ## invalid nrow and ncol
  new.cp <- lapply(old.cp, as, "CsparseMatrix")
  new.cp[[1]] <- new.cp[[1]][1:10, ]
  expect_error(colPairs(sce) <- new.cp, "when coerce to SelfHits")

  ## valid nrow and ncol but wrong nnode
  new.cp <- lapply(old.cp, as, "CsparseMatrix")
  new.cp[[1]] <- new.cp[[1]][1:10, 1:10]
  expect_error(colPairs(sce) <- new.cp)

  ## valid nnode but wrong names
  new.cp <- lapply(old.cp, as, "CsparseMatrix")
  colnames(new.cp[[1]]) <- rownames(new.cp[[1]]) <- seq_len(ncol(new.cp[[1]]))
  expect_error(colPairs(sce) <- new.cp)
})

test_that("colPair<- type = 'missing'", {
  sce <- exampleSCFE()

  old.cp <- colPairs(sce)

  ## New colPair is a SelfHits
  colPair(sce) <- old.cp[[2]]
  expect_equal(colPair(sce), old.cp[[2]])
  expect_equal(colPair(sce, 1), old.cp[[2]])
  expect_equal(colPair(sce, colPairNames(sce)[[1]]), old.cp[[2]])

  ## New colPair is a sparse matrix
  v <- as(old.cp[[2]], "CsparseMatrix")
  colPair(sce) <- v
  expect_equal(unname(colPair(sce, asSparse = TRUE)), v)
  expect_equal(unname(colPair(sce, 1, asSparse = TRUE)), v)
  expect_equal(unname(colPair(sce, colPairNames(sce)[[1]], asSparse = TRUE)), v)
})

test_that("colPair<- type = 'numeric'", {
  sce <- exampleSCFE()

  old.cp <- colPairs(sce)

  ## New colPair is a SelfHits
  colPair(sce, 2) <- old.cp[[1]]
  expect_equal(colPair(sce, 2), old.cp[[1]])
  expect_equal(colPair(sce, colPairNames(sce)[[2]]), old.cp[[1]])
  expect_equal(colPairNames(sce), names(old.cp))

  ## New colPair is a sparse matrix
  v <- as(old.cp[[1]], "CsparseMatrix")
  colPair(sce, 2) <- v
  expect_equal(unname(colPair(sce, 2, asSparse = TRUE)), v)
  expect_equal(unname(colPair(sce, colPairNames(sce)[[2]], asSparse = TRUE)), v)

  ## Cannot add a new colPair with numeric index
  expect_error(colPair(sce, 3) <- old.cp[[2]])

  ## Remove all colPairs first. Still cannot add a new colPair
  colPairs(sce) <- NULL
  expect_error(colPair(sce, 1) <- old.cp[[2]])
})

test_that("colPair<- type = 'character'", {
  sce <- exampleSCFE()

  old.cp <- colPairs(sce)

  colPair(sce, colPairNames(sce)[[1]]) <- old.cp[[2]]
  expect_equal(colPair(sce, 1), old.cp[[2]])
  expect_equal(colPair(sce, colPairNames(sce)[[1]]), old.cp[[2]])
  expect_equal(colPairNames(sce), names(old.cp))

  v <- as(old.cp[[2]], "CsparseMatrix")
  colPair(sce, colPairNames(sce)[[1]]) <- v
  expect_equal(unname(colPair(sce, 1, asSparse = TRUE)), v)
  expect_equal(unname(colPair(sce, colPairNames(sce)[[1]], asSparse = TRUE)), v)

  ## Can add a new colPair with character index
  colPair(sce, "new") <- old.cp[[2]]
  expect_equal(colPair(sce, 3), old.cp[[2]])
  expect_equal(colPair(sce, colPairNames(sce)[[3]]), old.cp[[2]])
  expect_equal(colPairNames(sce), c(names(old.cp), "new"))

  ## Remove all colPairs, then add a new colPair
  colPairs(sce) <- NULL
  colPair(sce, "new") <- old.cp[[2]]
  expect_equal(colPair(sce), old.cp[[2]])
  expect_equal(colPair(sce, 1), old.cp[[2]])
  expect_equal(colPair(sce, colPairNames(sce)[[1]]), old.cp[[2]])
  expect_equal(colPairNames(sce), "new")
})

test_that("colPairNames<-", {
  sce <- exampleSCFE()

  old.cp <- colPairs(sce)

  ## Modify an existing colPairNames
  colPairNames(sce)[2] <- "SNN2"
  expect_equal(colPairNames(sce)[2], "SNN2")
  expect_equal(colPairNames(sce)[1], names(old.cp)[1])

  ## Remove all colPairNames
  expect_warning(colPairNames(sce) <- character())
  expect_equal(colPairNames(sce), paste0("unnamed", seq(colPairNames(sce))))
  expect_error(colPairNames(sce) <- NULL)
})

test_that("rowPairs getters", {
  sce <- exampleSCFE()
  sce <- add_rowPair(sce)

  get_rowPair_test(sce)
  get_rowPair_test(sce, 1)
  get_rowPair_test(sce, 1L)
  get_rowPair_test(sce, "cor1")
  get_rowPair_test(sce, 2)
  get_rowPair_test(sce, 2L)
  get_rowPair_test(sce, "cor2")

  all.rp <- rowPairs(sce)
  expect_s4_class(all.rp, "SimpleList")
  expect_equal_no_attr(names(all.rp), rowPairNames(sce))
  expect_equal_no_attr(
    as.list(all.rp),
    lapply(rowPairNames(sce), rowPair, x = sce)
  )

  ## `asSparse = TRUE` works
  all.rp <- rowPairs(sce, asSparse = TRUE)
  expect_s4_class(all.rp, "SimpleList")
  expect_equal_no_attr(names(all.rp), rowPairNames(sce))
  expect_equal_no_attr(
    as.list(all.rp),
    lapply(rowPairNames(sce), rowPair, x = sce, asSparse = TRUE)
  )

  expect_equal_no_attr(rowPair(sce), rowPairs(sce)[[1]])
  expect_equal_no_attr(rowPair(sce, 2), rowPairs(sce)[["cor2"]])
  expect_equal_no_attr(rowPair(sce, "cor2"), rowPairs(sce)[[2]])

  expect_equal_no_attr(
    rowPair(sce, asSparse = TRUE),
    rowPairs(sce, asSparse = TRUE)[[1]]
  )
  expect_equal_no_attr(
    rowPair(sce, 2, asSparse = TRUE),
    rowPairs(sce, asSparse = TRUE)[["cor2"]]
  )
  expect_equal_no_attr(
    rowPair(sce, "cor2", asSparse = TRUE),
    rowPairs(sce, asSparse = TRUE)[[2]]
  )
})

test_that("rowPairs<- set NULL", {
  sce <- exampleSCFE()
  sce <- add_rowPair(sce)

  rowPairs(sce) <- NULL
  expect_in("rowPairs", colnames(int_elementMetadata(sce)))
  expect_length(rowPairs(sce), 0)
  expect_length(rowPairNames(sce), 0)
})

test_that("rowPairs<- set new list", {
  sce <- exampleSCFE()
  sce <- add_rowPair(sce)

  old.rp <- rowPairs(sce)
  rowPairs(sce) <- old.rp
  expect_equal_no_attr(rowPairs(sce), old.rp)

  ## new rowPairs without names
  new.rp <- unname(old.rp)
  expect_warning(rowPairs(sce) <- new.rp, "unnamed")
  expect_equal(rowPairNames(sce), paste0("unnamed", seq(new.rp)))

  new.rp <- list(cor1 = old.rp[["cor1"]], old.rp[[2]])
  expect_warning(rowPairs(sce) <- new.rp, "unnamed")
  expect_equal(rowPairNames(sce)[2], "unnamed2")
  expect_equal(rowPairNames(sce)[1], names(old.rp)[1])

  ## new rowPairs with wrong nnode
  new.rp <- old.rp
  new.rp[[1]] <- PairedHits(new.rp[[1]])[1:10]@hits
  expect_error(rowPairs(sce) <- new.rp)

  ## new rowPairs as sparseMatrix
  ## invalid nrow and ncol
  new.rp <- lapply(old.rp, as, "CsparseMatrix")
  new.rp[[1]] <- new.rp[[1]][1:10, ]
  expect_error(rowPairs(sce) <- new.rp, "when coerce to SelfHits")

  ## valid nrow and ncol but wrong nnode
  new.rp <- lapply(old.rp, as, "CsparseMatrix")
  new.rp[[1]] <- new.rp[[1]][1:10, 1:10]
  expect_error(rowPairs(sce) <- new.rp)

  ## valid nnode but wrong names
  new.rp <- lapply(old.rp, as, "CsparseMatrix")
  colnames(new.rp[[1]]) <- rownames(new.rp[[1]]) <- seq_len(ncol(new.rp[[1]]))
  expect_error(rowPairs(sce) <- new.rp)
})

test_that("rowPair<- type = 'missing'", {
  sce <- exampleSCFE()
  sce <- add_rowPair(sce)

  old.rp <- rowPairs(sce)

  ## New rowPair is a SelfHits
  rowPair(sce) <- old.rp[[2]]
  expect_equal(rowPair(sce), old.rp[[2]])
  expect_equal(rowPair(sce, 1), old.rp[[2]])
  expect_equal(rowPair(sce, rowPairNames(sce)[[1]]), old.rp[[2]])

  ## New rowPair is a sparse matrix
  v <- as(old.rp[[2]], "CsparseMatrix")
  rowPair(sce) <- v
  expect_equal(unname(rowPair(sce, asSparse = TRUE)), v)
  expect_equal(unname(rowPair(sce, 1, asSparse = TRUE)), v)
  expect_equal(unname(rowPair(sce, rowPairNames(sce)[[1]], asSparse = TRUE)), v)
})

test_that("rowPair<- type = 'numeric'", {
  sce <- exampleSCFE()
  sce <- add_rowPair(sce)

  old.rp <- rowPairs(sce)

  ## New rowPair is a SelfHits
  rowPair(sce, 2) <- old.rp[[1]]
  expect_equal(rowPair(sce, 2), old.rp[[1]])
  expect_equal(rowPair(sce, rowPairNames(sce)[[2]]), old.rp[[1]])
  expect_equal(rowPairNames(sce), names(old.rp))

  ## New rowPair is a sparse matrix
  v <- as(old.rp[[1]], "CsparseMatrix")
  rowPair(sce, 2) <- v
  expect_equal(unname(rowPair(sce, 2, asSparse = TRUE)), v)
  expect_equal(unname(rowPair(sce, rowPairNames(sce)[[2]], asSparse = TRUE)), v)

  ## Cannot add a new rowPair with numeric index
  expect_error(rowPair(sce, 3) <- old.rp[[2]])

  ## Remove all rowPairs first. Still cannot add a new rowPair
  rowPairs(sce) <- NULL
  expect_error(rowPair(sce, 1) <- old.rp[[2]])
})

test_that("rowPair<- type = 'character'", {
  sce <- exampleSCFE()
  sce <- add_rowPair(sce)

  old.rp <- rowPairs(sce)

  rowPair(sce, rowPairNames(sce)[[1]]) <- old.rp[[2]]
  expect_equal(rowPair(sce, 1), old.rp[[2]])
  expect_equal(rowPair(sce, rowPairNames(sce)[[1]]), old.rp[[2]])
  expect_equal(rowPairNames(sce), names(old.rp))

  v <- as(old.rp[[2]], "CsparseMatrix")
  rowPair(sce, rowPairNames(sce)[[1]]) <- v
  expect_equal(unname(rowPair(sce, 1, asSparse = TRUE)), v)
  expect_equal(unname(rowPair(sce, rowPairNames(sce)[[1]], asSparse = TRUE)), v)

  ## Can add a new rowPair with character index
  rowPair(sce, "new") <- old.rp[[2]]
  expect_equal(rowPair(sce, 3), old.rp[[2]])
  expect_equal(rowPair(sce, rowPairNames(sce)[[3]]), old.rp[[2]])
  expect_equal(rowPairNames(sce), c(names(old.rp), "new"))

  ## Remove all rowPairs, then add a new rowPair
  rowPairs(sce) <- NULL
  rowPair(sce, "new") <- old.rp[[2]]
  expect_equal(rowPair(sce), old.rp[[2]])
  expect_equal(rowPair(sce, 1), old.rp[[2]])
  expect_equal(rowPair(sce, rowPairNames(sce)[[1]]), old.rp[[2]])
  expect_equal(rowPairNames(sce), "new")
})

test_that("rowPairNames<-", {
  sce <- exampleSCFE()
  sce <- add_rowPair(sce)

  old.rp <- rowPairs(sce)

  ## Modify an existing rowPairNames
  rowPairNames(sce)[2] <- "cor2"
  expect_equal(rowPairNames(sce)[2], "cor2")
  expect_equal(rowPairNames(sce)[1], names(old.rp)[1])

  ## Remove all rowPairNames
  expect_warning(rowPairNames(sce) <- character())
  expect_equal(rowPairNames(sce), paste0("unnamed", seq(rowPairNames(sce))))
  expect_error(rowPairNames(sce) <- NULL)
})

test_that("SingleCellFlexExperiment subset", {
  pbmc <- exampleSCFE()

  # integer indices
  i <- sample(1:nrow(pbmc), 20)
  j <- sample(1:ncol(pbmc), 10)
  pbmc2 <- test_subset_SCE(pbmc, i = i, j = j)
  pbmc2 <- test_subset_SCE(pbmc, i = i)
  pbmc2 <- test_subset_SCE(pbmc, j = j)

  # character indices
  i <- sample(rownames(pbmc), 10)
  j <- sample(colnames(pbmc))
  pbmc2 <- test_subset_SCE(pbmc, i = i, j = j)
  pbmc2 <- test_subset_SCE(pbmc, i = i)
  pbmc2 <- test_subset_SCE(pbmc, j = j)

  # logical indices
  i <- rowData(pbmc)$vst.mean > 0.5
  j <- pbmc$nCount_RNA > 20
  pbmc2 <- test_subset_SCE(pbmc, i = i, j = j)
  pbmc2 <- test_subset_SCE(pbmc, i = i)
  pbmc2 <- test_subset_SCE(pbmc, j = j)
})






test_subset_RSE <- function(obj, i = NULL, j = NULL) {
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
  if (!is.null(i)) {
    expect_identical(rr[i], rr2)
    expect_identical(rd[i, , drop = FALSE], rd2)
  }
  if (!is.null(j)) {
    expect_identical(cd[j, , drop = FALSE], cd2)
  }
  return(obj2)
}

test_that("FlexExperiment .DollarNames", {

  fe <- exampleFE()

  ## .DollarNames should get all column names of colData(x)
  dollar.names <- .DollarNames(fe)
  expect_identical(dollar.names, colnames(colData(fe)))

  ## Use `$` to get or set one column
  new_ID <- paste0("AAA_", colnames(fe))
  fe$new_ID <- new_ID
  expect_identical(fe$new_ID, new_ID)

  nUMI <- Matrix::colSums(assay(fe))
  fe$nUMI <- nUMI
  expect_equal(fe$nUMI, nUMI)

  ## Use `[[]]` to get or set one column
  nGene <- Matrix::colSums(assay(fe) > 0)
  fe[["nGene"]] <- nGene
  expect_equal(fe[["nGene"]], nGene)
  expect_equal(fe$nGene, nGene)
})

test_that("FlexExperiment 'names'", {

  fe <- exampleFE()

  expect_equal_no_attr(names(fe), names(assays(fe)))

  ## Set new names
  fe2 <- fe
  new_names <- paste0("AAA_", names(fe))
  names(fe2) <- new_names
  expect_equal_no_attr(names(fe2), new_names)
  expect_equal_no_attr(names(assays(fe2)), new_names)

  ## Set empty names
  fe2 <- fe
  names(fe2) <- NULL
  expect_equal_no_attr(names(fe2), names(assays(fe2)))
  expect_null(names(fe2))
})

test_that("FlexExperiment 'rownames'", {

  fe <- exampleFE()

  expect_equal_no_attr(rownames(fe), rownames(fe@assays))
  expect_equal(dim(fe), dim(fe@assays))

  ## Set new rownames
  fe2 <- fe
  new_names <- paste0("AAA_", rownames(fe))
  rownames(fe2) <- new_names
  expect_equal_no_attr(rownames(fe2@assays), new_names)
  expect_equal_no_attr(rownames(fe2), new_names)

  ## Set empty rownames should fail for non-empty FlexExperiment
  fe2 <- fe
  expect_warning(rownames(fe2) <- NULL)
  expect_equal_no_attr(rownames(fe2), rownames(fe))
})

test_that("FlexExperiment 'colnames'", {

  fe <- exampleFE()

  expect_equal_no_attr(colnames(fe), colnames(fe@assays))
  expect_equal(dim(fe), dim(fe@assays))

  ## Set new colnames
  fe2 <- fe
  new_names <- paste0("AAA_", colnames(fe))
  colnames(fe2) <- new_names
  expect_equal_no_attr(colnames(fe2@assays), new_names)
  expect_equal_no_attr(colnames(fe2), new_names)

  ## Set empty colnames should fail for non-empty FlexExperiment
  fe2 <- fe
  expect_warning(colnames(fe2) <- NULL)
  expect_equal_no_attr(colnames(fe2), colnames(fe))
})

test_that("FlexExperiment 'dimnames'", {

  fe <- exampleFE()

  ## dimnames
  expect_equal_no_attr(dimnames(fe), list(rownames(fe), colnames(fe)))
  expect_equal_no_attr(dimnames(fe@assays), list(rownames(fe), colnames(fe)))

  ## Set new dimnames
  fe2 <- fe
  new_dimnames <- lapply(dimnames(fe), paste0, "_AAA")
  dimnames(fe2) <- new_dimnames
  expect_equal_no_attr(dimnames(fe2), new_dimnames)
  expect_equal_no_attr(dimnames(fe2@assays), new_dimnames)

  ## Set new dimnames with empty colnames
  fe2 <- fe
  new_names <- paste0(rownames(fe), "_AAA")
  expect_warning(dimnames(fe2) <- list(new_names, NULL))
  expect_equal_no_attr(dimnames(fe2), list(new_names, colnames(fe)))

  ## Set new dimnames with empty rownames
  fe2 <- fe
  new_names <- paste0(colnames(fe), "_AAA")
  expect_warning(dimnames(fe2) <- list(NULL, new_names))
  expect_equal_no_attr(dimnames(fe2), list(rownames(fe), new_names))

  ## Set empty dimnames should fail
  fe2 <- fe
  expect_warning(dimnames(fe2) <- NULL)
  expect_equal_no_attr(dimnames(fe2), dimnames(fe))
  expect_warning(dimnames(fe2) <- list())
  expect_equal_no_attr(dimnames(fe2), dimnames(fe))

  ## Set dimnames with wrong length
  expect_error(dimnames(fe2) <- list(rownames(fe)))
  expect_error(dimnames(fe2) <- list(colnames(fe)))
  expect_error(dimnames(fe2) <- list(rownames(fe), NULL, colnames(fe)))
})

test_that("FlexExperiment 'rowData'", {

  fe <- exampleFE()

  expect_equal_no_attr(rownames(rowData(fe)), rownames(fe))

  ## Set new `rownames(x)` will affect `rownames(rowData(x))`
  fe2 <- fe
  new_ID <- paste0("AAA_", rownames(fe2))
  rownames(fe2) <- new_ID
  expect_equal_no_attr(rownames(rowData(fe2)), new_ID)
  expect_equal_no_attr(rownames(fe2), new_ID)

  ## Set new `rowData` with re-ordered rows will force to match original
  ## `rownames(x)`
  fe2 <- fe
  rdata <- rowData(fe2)
  rowData(fe2) <- rdata[sample(rownames(rdata)), ]
  expect_identical(rowData(fe2), rdata)

  ## Set new `rowData` in base data.frame
  fe2 <- fe
  rdata <- rowData(fe2)
  rowData(fe2) <- as.data.frame(rdata[sample(rownames(rdata)), ])
  expect_identical(rowData(fe2), rdata)

  ## Remove all columns in `rowData(x)`
  fe2 <- fe
  rowData(fe2) <- NULL
  expect_equal(ncol(rowData(fe2)), 0)
  expect_equal(nrow(rowData(fe2)), nrow(fe2))

  ## Modify one column in `rowData(x)`
  fe2 <- fe
  nCells <- Matrix::rowSums(assay(fe2) > 0)
  rowData(fe2)[['nCells']] <- nCells
  expect_equal_no_attr(rowData(fe2)[['nCells']], nCells)

  ## Modify several columns in `rowData(x)`
  fe2 <- fe
  rowData(fe2)[, c("GeneID", "nCells")] <- DataFrame(
    GeneID = rowData(fe2)[['GeneName']],
    nCells = nCells
  )
  expect_equal_no_attr(rowData(fe2)[['nCells']], nCells)
  expect_equal_no_attr(rowData(fe2)[['GeneID']], rowData(fe2)[['GeneName']])

  ## Cannot set `rowData` without row names
  rdata <- rowData(fe, use.names = FALSE)
  expect_null(rownames(rdata))
  expect_error(rowData(fe) <- rdata)

  ## Cannot set `rowData` with unmatched row names
  fe2 <- fe
  rdata <- rowData(fe2)
  expect_error(rowData(fe2) <- rdata[sample(rownames(rdata), 10), ])

  rdata2 <- rdata
  rownames(rdata2)[1:10] <- letters[1:10]
  expect_error(rowData(fe2) <- rdata2)

  ## rowData(x, use.names = FALSE) works
  rdata <- rowData(fe, use.names = FALSE)
  expect_null(rownames(rdata))
  expect_equal(nrow(rdata), nrow(fe))
})

test_that("FlexExperiment 'rowRanges'", {

  fe <- exampleFE()

  ## names(rowRanges(x)) should be the same as rownames(x)
  new_ID <- paste0("xxxx", rownames(fe))
  rownames(fe) <- new_ID
  expect_equal_no_attr(rownames(rowData(fe)), new_ID)
  expect_equal_no_attr(names(rowRanges(fe)), new_ID)

  ## Set new `rowRanges` reordered will force to match original rownames(x)
  old_rr <- rowRanges(fe)
  rowRanges(fe) <- old_rr[sample(names(old_rr))]
  expect_identical(rowRanges(fe), old_rr)

  ## Cannot remove `rowRanges` since we need to keep `rowData`
  old_rr <- rowRanges(fe)
  expect_warning(rowRanges(fe) <- NULL)
  expect_identical(rowRanges(fe), old_rr)

  ## Cannot set new `rowRanges` without feature names
  new_rr <- old_rr
  names(new_rr) <- NULL
  expect_error(rowRanges(fe) <- new_rr)

  ## Cannot set new `rowRanges` with unmatched feature names
  new_rr <- old_rr
  expect_error(rowRanges(fe) <- new_rr[1:10])

  names(new_rr)[1:10] <- paste0("xxxx", names(old_rr)[1:10])
  expect_error(rowRanges(fe) <- new_rr)

  ## Set new `rowRanges` without metadata will remove original `rowData`
  new_rr <- old_rr
  mcols(new_rr) <- NULL
  rowRanges(fe) <- new_rr
  expect_equal(ncol(rowData(fe)), 0)
})

test_that("FlexExperiment 'colData'", {

  fe <- exampleFE()

  expect_equal_no_attr(rownames(colData(fe)), colnames(fe))

  ## Set new `colnames(x)` will affect `rownames(colData(x))`
  new_ID <- paste0("AAA_", colnames(fe))
  colnames(fe) <- new_ID
  expect_equal_no_attr(rownames(colData(fe)), new_ID)
  expect_equal_no_attr(colnames(fe), new_ID)

  ## Set new `colData` with re-ordered rows will force to match original
  ## `colnames(x)`
  fe2 <- fe
  cdata <- colData(fe2)
  colData(fe2) <- cdata[sample(rownames(cdata)), ]
  expect_identical(colData(fe2), cdata)

  ## Modify one column in `colData(x)`
  fe2 <- fe
  nUMI <- Matrix::colSums(assay(fe2))
  nGene <- Matrix::colSums(assay(fe2) > 0)
  colData(fe2)[['nGene']] <- nGene
  expect_equal_no_attr(colData(fe2)[['nGene']], nGene)

  ## Modify several columns in `rowData(x)`
  fe2 <- fe
  colData(fe2)[, c("nGene", "nUMI")] <- DataFrame(nGene, nUMI)
  expect_equal_no_attr(colData(fe2)[['nGene']], nGene)
  expect_equal_no_attr(colData(fe2)[['nUMI']], nUMI)

  ## Set new `colData` in base data.frame
  cdata <- colData(fe)
  colData(fe) <- as.data.frame(cdata[sample(rownames(cdata)), ])
  expect_identical(colData(fe), cdata)

  ## Remove all columns in `colData(x)`
  fe2 <- fe
  colData(fe2) <- NULL
  expect_equal(ncol(colData(fe2)), 0)
  expect_equal(nrow(colData(fe2)), ncol(fe2))

  ## Cannot set `colData` without row names
  cdata <- colData(fe)
  rownames(cdata) <- NULL
  expect_error(cowData(fe) <- cdata)

  ## Cannot set `colData` with unmatched row names
  fe2 <- fe
  cdata <- colData(fe2)
  expect_error(colData(fe2) <- cdata[sample(rownames(cdata), 10), ])

  cdata2 <- cdata
  rownames(cdata2)[1:10] <- letters[1:10]
  expect_error(colData(fe2) <- cdata2)
})

test_that("FlexExperiment 'assays'", {

  fe <- exampleFE()

  ## assays
  expect_s4_class(assays(fe), "SimpleList")
  expect_length(fe, length(assays(fe)))
  expect_length(fe@assays, length(assays(fe)))
  expect_identical(names(fe), names(assays(fe)))
  expect_identical(names(fe@assays), names(assays(fe)))

  mats <- assays(fe)
  expect_equal(mats[[1]], assay(fe))
  expect_equal(mats[[2]], assay(fe, 2))
  expect_equal(mats[[2]], assay(fe, "logcounts"))
})

test_that("FlexExperiment 'assays<-'", {

  fe <- exampleFE()

  ## Cannot set `assays` with empty dimension names
  fe2 <- fe
  new.assays <- assays(fe2, withDimnames = FALSE)
  expect_error(assays(fe2) <- new.assays)

  ## Set `assays` with one assay containing empty dimension names
  new.assays <- assays(fe2)
  dimnames(new.assays[[1]]) <- NULL
  assays(fe2) <- new.assays
  expect_equal(dim(assay(fe2)), c(0, 0))

  ## Some new assays are subset
  ## The new x@assays must match the old rowData and colData
  fe2 <- fe
  new.assays <- assays(fe2)
  new.assays[[1]] <- new.assays[[1]][1:10, 1:20]
  new.assays[[2]] <- new.assays[[2]][11:30, 31:65]
  new.assays[[3]] <- new.assays[[3]][91:100, 61:70]
  expect_error(assays(fe2) <- new.assays)

  ## `withDimnames = FALSE` should not work
  fe2 <- fe
  new.assays <- assays(fe2)
  dimnames(new.assays[[1]]) <- NULL
  expect_warning(assays(fe2, withDimnames = FALSE) <- new.assays)
  expect_equal(dim(assay(fe2)), c(0, 0))

  fe2 <- fe
  expect_warning(assay(fe2, withDimnames = FALSE) <- assay(fe2)[1:10, 1:10])
  expect_equal(assay(fe2), assay(fe)[1:10, 1:10])
})

test_that("FlexExperiment subset", {

  fe <- exampleFE()

  # integer indices
  i <- sample(1:nrow(fe), 10)
  j <- sample(1:ncol(fe), 20)
  fe2 <- test_subset_RSE(fe, i = i, j = j)
  fe2 <- test_subset_RSE(fe, i = i)
  fe2 <- test_subset_RSE(fe, j = j)

  # character indices
  i <- sample(rownames(fe), 10)
  j <- sample(colnames(fe))
  fe2 <- test_subset_RSE(fe, i = i, j = j)
  fe2 <- test_subset_RSE(fe, i = i)
  fe2 <- test_subset_RSE(fe, j = j)

  # logical indices
  i <- rowData(fe)$nCells > 3
  j <- fe$nCount > 20
  fe2 <- test_subset_RSE(fe, i = i, j = j)
  fe2 <- test_subset_RSE(fe, i = i)
  fe2 <- test_subset_RSE(fe, j = j)
})

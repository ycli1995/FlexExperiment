
generate_random_genes <- function(n = 10, max = round(1.2 * n)) {
  number <- sample(seq_len(max), n)
  paste0("gene", number)
}

generate_random_barcodes <- function(n = 10, length = 8) {
  bases <- c("A", "C", "G", "T")
  replicate(n, paste0(sample(bases, length, replace = TRUE), collapse = ""))
}

generate_sparse_matrix <- function(nrow, ncol, density = 0.75) {
  mat1 <- Matrix::rsparsematrix(nrow = nrow, ncol = ncol, density = density)
  rownames(mat1) <- generate_random_genes(nrow(mat1))
  colnames(mat1) <- generate_random_barcodes(ncol(mat1))
  mat1
}

test_concat_both <- function(...) {
  mat.list <- list(...)

  # concatByRows
  m <- concatByRows(mat.list[[1]], mat.list[2:length(mat.list)])
  expect_equal(nrow(m), sum(sapply(mat.list, nrow)))
  expect_setequal(colnames(m), unique(unlist(lapply(mat.list, colnames))))
  expect_equal(anyDuplicated(rownames(m)), 0)
  expect_equal(anyDuplicated(colnames(m)), 0)

  offset <- 0
  for (i in seq_along(mat.list)) {
    expect_equal(
      unname(m[1:nrow(mat.list[[i]]) + offset, colnames(mat.list[[i]])]),
      unname(mat.list[[i]])
    )
    offset <- offset + nrow(mat.list[[i]])
    # Prepare for concatByCols
    mat.list[[i]] <- Matrix::t(mat.list[[i]])
  }

  # concatByCols
  m <- concatByCols(mat.list[[1]], mat.list[2:length(mat.list)])
  expect_equal(ncol(m), sum(sapply(mat.list, ncol)))
  expect_setequal(rownames(m), unique(unlist(lapply(mat.list, rownames))))
  expect_equal(anyDuplicated(rownames(m)), 0)
  expect_equal(anyDuplicated(colnames(m)), 0)

  offset <- 0
  for (i in seq_along(mat.list)) {
    expect_equal(
      unname(m[rownames(mat.list[[i]]), 1:ncol(mat.list[[i]]) + offset]),
      unname(mat.list[[i]])
    )
    offset <- offset + ncol(mat.list[[i]])
  }
  invisible(NULL)
}

test_that("concatenation for matrices", {
  m1 <- as.matrix(generate_sparse_matrix(8, 6))
  m2 <- as.matrix(generate_sparse_matrix(6, 9))
  m3 <- as.matrix(generate_sparse_matrix(8, 6))

  test_concat_both(m1, m2, m3)

  # duplicated row names
  rownames(m3)[1:4] <- rownames(m1)[5:8]
  test_concat_both(m1, m2, m3)

  # duplicated column names will raise an error for concatByRows
  colnames(m3)[4] <- colnames(m3)[1]
  expect_error(concatByRows(m1, list(m2, m3)))

  # duplicated row names will raise an error for concatByCols
  rownames(m3)[4] <- rownames(m3)[1]
  expect_error(concatByCols(m1, list(m2, m3)))
})

test_that("concatenation for sparse matrices", {
  m1 <- generate_sparse_matrix(8, 6)
  m2 <- generate_sparse_matrix(6, 9)
  m3 <- generate_sparse_matrix(8, 6)

  test_concat_both(m1, m2, m3)

  # duplicated row names
  rownames(m3)[1:4] <- rownames(m1)[5:8]
  test_concat_both(m1, m2, m3)

  # duplicated column names will raise an error for concatByRows
  colnames(m3)[4] <- colnames(m3)[1]
  expect_error(concatByRows(m1, list(m2, m3)))

  # duplicated row names will raise an error for concatByCols
  rownames(m3)[4] <- rownames(m3)[1]
  expect_error(concatByCols(m1, list(m2, m3)))
})


generate_NN <- function(n_obs, k) {
  idx <- matrix(0L, n_obs, k)
  for (i in seq_len(nrow(idx))) {
    idx[i, ] <- sample(seq_len(n_obs), k)
  }
  dist <- matrix(runif(n_obs * k), n_obs, k)
  rownames(idx) <- rownames(dist) <- as.character(seq_len(nrow(idx)))
  return(list(idx = idx, dist = dist))
}

generate_PairedHits <- function(n_obs, k) {
  nn <- generate_NN(n_obs, k)
  nnToPairedHits(nn$idx, nn$dist)
}

test_that("PairedHits", {
  ph <- generate_PairedHits(10, 5)

  expect_equal(length(ph), 10)
  expect_identical(names(ph), as.character(seq_along(ph)))

  names(ph) <- paste0("obs", names(ph))
  expect_identical(names(ph), paste0("obs", seq_len(length(ph))))

  # Subset a PairedHits by logical
  idx1 <- sample(c(TRUE, FALSE), length(ph), replace = TRUE)
  nm1 <- names(ph)[idx1]
  ph1 <- ph[idx1]
  expect_equal(names(ph1), nm1)

  # Subset a PairedHits by names
  nm1 <- sample(names(ph), length(5))
  ph1 <- ph[nm1]
  expect_equal(names(ph1), nm1)

  # Subset a PairedHits by indices
  idx1 <- sample(seq_len(length(ph)), length(5))
  nm1 <- names(ph)[idx1]
  ph1 <- ph[idx1]
  expect_equal(names(ph1), nm1)
})

test_that("PairedHits to sparse matrix", {
  ph <- generate_PairedHits(10, 5)
  ph2 <- as(ph, "CsparseMatrix")
  ph3 <- as(as(ph2, "PairedHits"), "CsparseMatrix")
  expect_identical(ph2, ph3)
})

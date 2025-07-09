
test_that("LinearEmbeddingMatrix to matrix", {
  rd1 <- reducedDim(exampleSCE())
  mat1 <- as(rd1, "matrix")
  rd2 <- as(mat1, "LinearEmbeddingMatrix")

  expect_identical(as(rd2, "matrix"), mat1)
  expect_identical(colnames(rd2), colnames(mat1))
  expect_identical(rownames(rd2), rownames(mat1))
})

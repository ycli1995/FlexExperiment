#' @importFrom S4Vectors new2
#' @importFrom methods setAs
#'
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
NULL

setAs("matrix", "LinearEmbeddingMatrix", function(from) {
  sample.names <- rownames(from)
  factor.names <- colnames(from)
  from <- unname(from)
  loadings <- matrix(0, ncol = ncol(from), nrow = 0)
  factor.data <- new2("DFrame", nrows = ncol(from), check = FALSE)
  rownames(factor.data) <- factor.names
  lem <- LinearEmbeddingMatrix(
    sampleFactors = from,
    featureLoadings = loadings,
    factorData = factor.data
  )
  rownames(lem) <- sample.names
  lem
})

#' @importFrom Matrix drop0 nnzero sparseMatrix
#' @importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importClassesFrom Matrix CsparseMatrix generalMatrix TsparseMatrix
NULL

#' Convert Between K Nearest Neighbors (KNN) and Graph
#'
#' Convert K Nearest Neighbors (KNN) from / to sparse graph.
#'
#' @param ... `r .dot_param`
#'
#' @name nn-coerce
NULL

#' @param nn.idx A \code{\link{matrix}} for observations x K neighbor indices.
#' @param nn.weights `NULL` (default) or a \code{\link{matrix}} for distances
#' corresponding to `nn.idx`.
#' @param repr Which format of the returned sparse matrix:
#' \itemize{
#' \item `"C"` for `r .doc_links("CsparseMatrix")`
#' \item `"T"` for `r .doc_links("TsparseMatrix")`
#' }
#' @param self.loops Whether to allow self-loops in the output sparse matrix.
#'
#' @returns
#' \itemize{
#' \item `nnToGraph`: Returns a sparse matrix indicating the KNN graph. Each
#' column must contain the same number (k) of non-zero values.
#' }
#'
#' @export
#' @rdname nn-coerce
nnToGraph <- function(
    nn.idx,
    nn.weights = NULL,
    repr = c("C", "T"),
    self.loops = TRUE,
    ...
) {
  repr <- match.arg(repr)
  mat <- .nn_to_triplets(nn.idx, nn.weights, self.loops = self.loops)
  mat <- sparseMatrix(
    i = mat$i,
    j = mat$j,
    x = mat$x,
    repr = repr,
    index1 = FALSE,
    ...
  )
  colnames(mat) <- rownames(nn.idx)
  rownames(mat) <- rownames(nn.idx)
  mat
}

#' @returns
#' \itemize{
#' \item `nnToPairedHits`: Returns a \code{\link{PairedHits}} indicating the KNN
#' graph. Each left node should contain k right nodes.
#' }
#'
#' @export
#' @rdname nn-coerce
nnToPairedHits <- function(nn.idx, nn.weights = NULL, self.loops = TRUE) {
  mat <- .nn_to_triplets(nn.idx, nn.weights, self.loops = self.loops)
  hits <- SelfHits(mat$i + 1L, mat$j + 1L, nnode = nrow(nn.idx))
  if (!is.null(nn.weights)) {
    mcols(hits)[["x"]] <- mat$x
  }
  PairedHits(hits, NAMES = rownames(nn.idx))
}

.nn_to_triplets <- function(nn.idx, nn.weights = NULL, self.loops = TRUE) {
  if (is.null(nn.weights)) {
    return(rcpp_knn_to_sparse_idx(nn_idx = nn.idx, self_loops = self.loops))
  }
  rcpp_knn_to_sparse_weights(
    nn_idx = nn.idx,
    nn_dist = nn.weights,
    self_loops = self.loops
  )
}

#' @param object An object to be converted into KNN indices and weights.
#'
#' @returns
#' \itemize{
#' \item `graphToNN` and `hitsToNN`: Returns a list containing two
#' \code{\link{matrix}}. `idx` indicates the KNN indices and `dist` indicates
#' the weights (distances) between observations and their nearest neighbors.
#' }
#'
#' @rdname nn-coerce
#' @export graphToNN
graphToNN <- function(object, ...) {
  UseMethod("graphToNN", object)
}

#' @rdname nn-coerce
#' @export
#' @method graphToNN sparseMatrix
graphToNN.sparseMatrix <- function(object, ...) {
  x <- as(x, "CsparseMatrix")
  if (!inherits(object, "generalMatrix")) {
    object <- as(object, "generalMatrix")
  }
  k <- unique(diff(object@p))
  if (length(k) > 1) {
    stop("Not all observations have an equal number of neighbors.")
  }
  nn.idx <- matrix(object@i + 1L, ncol = k, byrow = TRUE)
  nn.dist <- matrix(object@x, ncol = k, byrow = TRUE)
  rownames(nn.idx) <- rownames(nn.dist) <- colnames(object)
  return(list(idx = nn.idx, dist = nn.dist))
}

#' @rdname nn-coerce
#' @export
hitsToNN <- function(object, ...) {
  object <- hitsToMat(object, repr = "C", index1 = TRUE)
  graphToNN(object, ...)
}

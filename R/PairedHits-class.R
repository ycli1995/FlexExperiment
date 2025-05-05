#' @include utils.R
#'
#' @importFrom FlexAssays getErrors
#' @importFrom S4Vectors from to Hits nLnode normalizeSingleBracketSubscript
#' nRnode nnode SelfHits
NULL

#' The PairedHits Class
#'
#' A class for storing paired hit objects. This structure wraps a
#' `r .doc_links("SelfHits")` object and an optional character vector for node
#' names.
#'
#' @slot hits A `r .doc_links("SelfHits")` object.
#' @slot NAMES A character vector or `NULL`, containing node names for `hits`.
#'
#' @seealso `r .doc_links("SelfHits")`
#'
#' @name PairedHits
#' @docType class
#' @exportClass PairedHits
setClass("PairedHits", contains = "Annotated", slots = c(
  hits = "SelfHits",
  NAMES = "character_OR_NULL"
))

#' @param hits A `r .doc_links("SelfHits")` object.
#' @param NAMES A character vector or `NULL`, containing node names for `hits`.
#' The number of elements must be equal to `nnode(hits)`. Duplicated node names
#' are not allowed.
#'
#' @export
#' @rdname PairedHits
PairedHits <- function(hits, NAMES = NULL) {
  new("PairedHits", hits = sort(hits), NAMES = NAMES)
}

## validation ##################################################################

setValidity2("PairedHits", function(x) validPairedHits(x, immediate. = FALSE))

#' Validation for PairedHits
#'
#' Check whether a \code{\link{PairedHits}} is valid.
#'
#' @param x A \code{\link{PairedHits}} to check.
#' @param immediate. Logical indicating if the errors should be raised
#' immediately. Set `immediate. = FALSE` can be useful for context of
#' `r .doc_links("validObject")`.
#'
#' @seealso [getErrors()]
#'
#' @returns
#' Returns an invisible `NULL` if `x` is a valid \code{\link{PairedHits}}. If
#' `immediate. = TRUE`, this function will raise an error when `x` is invalid.
#' Otherwise, it will return a character vector that collects all errors.
#'
#' @export
validPairedHits <- function(x, immediate. = TRUE) {
  .valid_PairedHits_names(x, immediate. = immediate.)
}

.valid_PairedHits_names <- function(x, immediate. = TRUE) {
  if (length(names(x)) == 0) {
    return(invisible(NULL))
  }
  if (anyDuplicated(names(x))) {
    e <- sprintf("Duplicated 'names(<%s>)' are not allowed.", class(x)[1])
    return(getErrors(e, err, immediate.))
  }
  if (length(names(x)) == nLnode(x)) {
    return(invisible(NULL))
  }
  fmt <- "Length of 'names(<%s>)' (%i) differs to 'nLnode(<%s>)' (%i)"
  err <- sprintf(fmt, class(x)[1], length(names(x)), class(x)[1], nLnode(x))
  getErrors(err, immediate. = immediate.)
}

#' Methods for PairedHits
#'
#' A collection of accessor methods for retrieving or modifying basic
#' information from a \code{\link{PairedHits}} object.
#'
#' @param x,object A \code{\link{FlexAssays}} object
#' @param ... `r .dot_param`
#' @param value `r .val_param` Typically a character vector for node names.
#'
#' @name PairedHits-methods
NULL

#' @export
#' @rdname PairedHits-methods
showPairedHits <- function(object) {
  show_S4_title(object)
  coolcat("Node names(%d): %s\n", names(object))
  show(object@hits)
  invisible(NULL)
}

#' @export
#' @rdname PairedHits-methods
setMethod("show", "PairedHits", showPairedHits)

#' @returns
#' \itemize{
#' \item `length`: An integer scalar indicating the number of nodes.
#' }
#'
#' @export
#' @rdname PairedHits-methods
setMethod("length", "PairedHits", function(x) nLnode(x@hits))

#' @returns
#' \itemize{
#' \item `names`: A character vector specifying the node names.
#' }
#'
#' @export
#' @rdname PairedHits-methods
setMethod("names", "PairedHits", function(x) x@NAMES)

#' @returns
#' \itemize{
#' \item `names<-`: Returns `x` with updated node names.
#' }
#'
#' @export
#' @rdname PairedHits-methods
setMethod("names<-", c("PairedHits", "character"), function(x, value) {
  if (length(value) == 0) {
    x@NAMES <- NULL
    return(x)
  }
  x@NAMES <- value
  .valid_PairedHits_names(x)
  x
})

#' @export
#' @rdname PairedHits-methods
setMethod("names<-", c("PairedHits", "NULL"), function(x, value) {
  x@NAMES <- value
  x
})

#' @export
#' @rdname PairedHits-methods
setMethod("nLnode", "PairedHits", function(x) nLnode(x@hits))

#' @export
#' @rdname PairedHits-methods
setMethod("nRnode", "PairedHits", function(x) nRnode(x@hits))

#' @export
#' @rdname PairedHits-methods
setMethod("nnode", "PairedHits", function(x) nLnode(x@hits))

#' @param i,j Indices specifying elements to extract. `j` will be ignored, since
#' subsetting a \code{\link{PairedHits}} object with `i` will simultaneously
#' work on left and right nodes of the underlying `r .doc_links("SelfHits")`.
#' @param drop Ignored by subsetting. This is only used for matching the S4
#' generic function.
#'
#' @returns
#' \itemize{
#' \item `[`: A subset \code{\link{PairedHits}} object.
#' }
#'
#' @export
#' @rdname PairedHits-methods
setMethod("[", "PairedHits", function(x, i, j, ..., drop = FALSE) {
  if (missing(i)) {
    return(x)
  }
  if (anyNA(i) && is.numeric(i)) {
    i <- as.integer(i)
  } else {
    i <- normalizeSingleBracketSubscript(i, x)
  }
  NAMES <- x@NAMES[i]
  hits <- x@hits
  left <- match(from(hits), i, nomatch = 0L)
  right <- match(to(hits), i, nomatch = 0L)

  idx <- which(left > 0 & right > 0)
  hits2 <- SelfHits(left[idx], right[idx], nnode = length(i))
  mcols(hits2) <- mcols(hits)[idx, , drop = FALSE]

  initialize(x, hits = sort(hits2), NAMES = NAMES, metadata = metadata(x))
})

## coerce ######################################################################

#' Helpers for Coercing PairedHits
#'
#' Helper functions for coercing between \code{\link{PairedHits}} object and
#' sparse matrix.
#'
#' @param x An object to be coerced.
#' @param ... `r .dot_param`
#'
#' @name PairedHits-coerce
NULL

#' @returns
#' \itemize{
#' \item `hitsToMat`: Returns a sparse matrix where row indices are `from(x)`
#' and column indices are `to(x)`.
#' }
#'
#' @rdname PairedHits-coerce
#' @export hitsToMat
hitsToMat <- function(x, ...) UseMethod("hitsToMat", x)

#' @rdname PairedHits-coerce
#' @export
#' @method hitsToMat Hits
hitsToMat.Hits <- function(x, ...) {
  m <- mcols(x)
  if (!is.null(m) && ncol(m)) {
    val <- m[, 1]
  } else {
    fmt<- "No values found in 'mcols(<%s>)', filling matrix with TRUE instead."
    warning(sprintf(fmt, class(x)[1]), immediate. = TRUE, call. = FALSE)
    val <- TRUE
  }
  sparseMatrix(from(x), to(x), x = val, dims = c(nLnode(x), nRnode(x)), ...)
}

#' @rdname PairedHits-coerce
#' @export
#' @method hitsToMat PairedHits
hitsToMat.PairedHits <- function(x, ...) {
  out <- hitsToMat(x@hits, ...)
  colnames(out) <- rownames(out) <- names(x)
  out
}

#' @returns
#' \itemize{
#' \item `matToHits`: Returns a `r .doc_links("Hits")` object where left nodes
#' indicate row indices and right nodes indicate column indices.
#' }
#'
#' @rdname PairedHits-coerce
#' @export matToHits
matToHits <- function(x, ...) UseMethod("matToHits", x)

#' @rdname PairedHits-coerce
#' @export
#' @method matToHits default
matToHits.default <- function(x, self = TRUE, ...) {
  if (self) {
    if (nrow(x) != ncol(x)) {
      fmt <- "'nrow' (%i) must be equal to 'ncol' (%i) when coerce to SelfHits."
      stop(sprintf(fmt, nrow(x), ncol(x)))
    }
  }
  x <- as(x, "TsparseMatrix")
  if (!inherits(x, "generalMatrix")) {
    x <- as(x, "generalMatrix")
  }
  if (self) {
    hits <- SelfHits(x@i + 1L, x@j + 1L, nnode = nrow(x), ...)
  } else {
    hits <- Hits(x@i + 1L, x@j + 1L, nLnode = nrow(x), nRnode = ncol(x), ...)
  }
  if (is.logical(x@x)) {
    return(hits)
  }
  mcols(hits)[["x"]] <- x@x
  hits
}

#' @returns
#' \itemize{
#' \item `matToPairedHits`: Returns a \code{\link{PairedHits}} object. The input
#' `x` must contain identical row names and column names.
#' }
#'
#' @rdname PairedHits-coerce
#' @export matToPairedHits
matToPairedHits <- function(x, ...) UseMethod("matToPairedHits", x)

#' @rdname PairedHits-coerce
#' @export matToPairedHits
#' @method matToPairedHits default
matToPairedHits.default <- function(x, ...) {
  if (!.identical_fmatch(rownames(x), colnames(x))) {
    stop("rownames(x) must be identical to colnames(x) for 'matToPairedHits'.")
  }
  hits <- matToHits(x)
  PairedHits(hits, NAMES = rownames(x))
}

### setAs ######################################################################

setAs("SelfHits", "PairedHits", function(from) PairedHits(from))
setAs("PairedHits", "SelfHits", function(from) from@hits)

setAs("PairedHits", "CsparseMatrix", function(from) hitsToMat(from, repr = "C"))
setAs("PairedHits", "TsparseMatrix", function(from) hitsToMat(from, repr = "T"))
setAs("PairedHits", "RsparseMatrix", function(from) hitsToMat(from, repr = "R"))

setAs("Hits", "CsparseMatrix", function(from) hitsToMat(from, repr = "C"))
setAs("Hits", "TsparseMatrix", function(from) hitsToMat(from, repr = "T"))
setAs("Hits", "RsparseMatrix", function(from) hitsToMat(from, repr = "R"))

setAs("TsparseMatrix", "SelfHits", function(from) matToHits(from, self = TRUE))
setAs("CsparseMatrix", "SelfHits", function(from) matToHits(from, self = TRUE))
setAs("RsparseMatrix", "SelfHits", function(from) matToHits(from, self = TRUE))

setAs("TsparseMatrix", "Hits", function(from) matToHits(from))
setAs("CsparseMatrix", "Hits", function(from) matToHits(from))
setAs("RsparseMatrix", "Hits", function(from) matToHits(from))

setAs("TsparseMatrix", "PairedHits", function(from) matToPairedHits(from))
setAs("CsparseMatrix", "PairedHits", function(from) matToPairedHits(from))
setAs("RsparseMatrix", "PairedHits", function(from) matToPairedHits(from))

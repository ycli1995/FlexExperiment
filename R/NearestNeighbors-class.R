#' @importFrom S4Vectors coolcat parallel_slot_names setValidity2
#' @importFrom Matrix as.matrix
#' @importFrom methods setClass setMethod show
#' @importClassesFrom S4Vectors character_OR_NULL
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class ########################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The NearestNeighbors class
#'
#' The `NearestNeighbors` class is designed to store the results of nearest
#' neighbor searching algorithms.
#'
#' @name NearestNeighbors
#' @docType class
#' @exportClass NearestNeighbors
setClass(
  "NearestNeighbors",
  contains = "Vector",
  slots = c(
    nnidx = "matrix",
    nndist = "matrix",
    alg.idx = "ANY",
    alg.info = "list",
    NAMES = "character_OR_NULL"
  )
)

## Constructor #################################################################

#' @param nnidx A \code{\link{matrix}} containing the nearest neighbor indices.
#' @param nndist A \code{\link{matrix}} containing the distances between those
#' observations in `nnidx`.
#' @param alg.idx A pre-built neighbor searching index. E.g. the `annoy` index.
#' @param alg.info Any information associated with the algorithm that may be
#' needed downstream.
#' @param NAMES Names of observations for which the nearest neighbors have been
#' computed.
#'
#' @examples
#' # Create a NearestNeighbors
#' nnidx <- matrix(ceiling(runif(100) * 10), nrow = 20, ncol = 5)
#' nndist <- matrix(runif(100), nrow = 20, ncol = 5)
#' nn <- NearestNeighbors(nnidx, nndist)
#'
#' @rdname NearestNeighbors
#' @export
NearestNeighbors <- function(
    nnidx,
    nndist = NULL,
    alg.idx = NULL,
    alg.info = list(),
    NAMES = NULL
) {
  NAMES <- NAMES %||% rownames(nnidx)
  rownames(nnidx) <- NULL
  nndist <- nndist %||% matrix(0, nrow = nrow(nnidx), ncol = 0)
  rownames(nndist) <- NULL
  new(
    "NearestNeighbors",
    nnidx = nnidx,
    nndist = nndist,
    alg.idx = alg.idx,
    alg.info = alg.info,
    NAMES = NAMES
  )
}

## Validity ####################################################################

validNearestNeighbors <- function(x, immediate. = TRUE) {
  err <- NULL
  if (!all(dim(x@nnidx) == dim(x@nndist))) {
    e <- sprintf(
      "`dim(nnidx)` (%s) must be equal to `dim(nndist)` (%s).",
      paste(dim(x@nnidx), collapse = ", "),
      paste(dim(x@nndist), collapse = ", ")
    )
    err <- getErrors(e, immediate. = immediate.)
  }
  if (length(x@NAMES) == 0) {
    return(err)
  }
  if (nrow(x@nnidx) != length(x@NAMES)) {
    e <- sprintf(
      "`length(NAMES)` (%d) must be equal to `nrow(nnidx)` (%d)",
      length(x@NAMES),
      nrow(x@nnidx)
    )
    err <- getErrors(e, err, immediate. = immediate.)
  }
  if (anyDuplicated(x@NAMES)) {
    e <- "Duplicated NAMES are not allowed."
    err <- getErrors(e, err, immediate. = immediate.)
  }
  invisible(err)
}

setValidity2("NearestNeighbors", function(x) validNearestNeighbors(x, FALSE))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NearestNeighbors-methods #####################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' NearestNeighbors methods
#'
#' Methods to manipulate a \code{\link{NearestNeighbors}} object
#'
#' @param object,x A \code{\link{NearestNeighbors}} object
#' @param ... `r .dot_param`
#'
#' @examples
#' # Create a NearestNeighbors
#' nnidx <- matrix(ceiling(runif(100) * 10), nrow = 20, ncol = 5)
#' nndist <- matrix(runif(100), nrow = 20, ncol = 5)
#' nn <- NearestNeighbors(nnidx, nndist)
#'
#' names(nn) <- letters[1:nrow(nn)]
#'
#' @name NearestNeighbors-methods
NULL

## Getters #####################################################################

#' @examples
#' nnidx(nn)
#'
#' @export
#' @rdname NearestNeighbors-methods
setGeneric("nnidx", function(x, ...) standardGeneric("nnidx"))

#' @param withDimnames Whether or not to add the names of observations.
#'
#' @export
#' @rdname NearestNeighbors-methods
setMethod("nnidx", "NearestNeighbors", function(x, withDimnames = TRUE, ...) {
  if (!withDimnames) {
    return(x@nnidx)
  }
  nnidx <- x@nnidx
  rownames(nnidx) <- x@NAMES
  nnidx
})

#' @examples
#' nndist(nn)
#'
#' @export
#' @rdname NearestNeighbors-methods
setGeneric("nndist", function(x, ...) standardGeneric("nndist"))

#' @param withDimnames Whether or not to add the names of observations.
#'
#' @export
#' @rdname NearestNeighbors-methods
setMethod("nndist", "NearestNeighbors", function(x, withDimnames = TRUE, ...) {
  if (!withDimnames) {
    return(x@nndist)
  }
  nndist <- x@nndist
  rownames(nndist) <- x@NAMES
  nndist
})


#' @export alg.idx
#' @rdname NearestNeighbors-methods
alg.idx <- function(object, ...) {
  UseMethod(generic = "alg.idx", object = object)
}

#' @export alg.info
#' @rdname NearestNeighbors-methods
alg.info <- function(object, ...) {
  UseMethod(generic = "alg.info", object = object)
}

#' @method alg.idx NearestNeighbors
#' @export
#' @rdname NearestNeighbors-methods
alg.idx.NearestNeighbors <- function(object, ...) {
  object@alg.idx
}

#' @method alg.info NearestNeighbors
#' @export
#' @rdname NearestNeighbors-methods
alg.info.NearestNeighbors <- function(object, ...) {
  object@alg.info
}

## as.matrix ###################################################################

#' @param use.dist Whether or not to use the distances for coercing into matrix.
#'
#' @examples
#' as.matrix(nn)
#' as.matrix(nn, use.dist = TRUE)
#'
#' @export
#' @rdname NearestNeighbors-methods
setMethod("as.matrix", "NearestNeighbors", function(x, use.dist = FALSE, ...) {
  if (use.dist) {
    return(nndist(x, ...))
  }
  nnidx(x, ...)
})

## show ########################################################################

#' @export
#' @rdname NearestNeighbors-methods
setMethod(
  f = "show",
  signature = "NearestNeighbors",
  definition = function(object) {
    cat("A", nrow(object), "x", ncol(object), class(object)[1], "\n")
    if (nrow(object) > 0) {
      coolcat("NAMES(%d): %s\n", names(object))
    }
    invisible(NULL)
  }
)

## names #######################################################################

#' @returns
#' \itemize{
#' \item `names`: A character vector specifying names of the observations.
#' }
#'
#' @examples
#' names(nn)
#'
#' @export
#' @rdname NearestNeighbors-methods
setMethod(
  f = "names",
  signature = "NearestNeighbors",
  definition = function(x) x@NAMES
)

#' @returns
#' \itemize{
#' \item `names<-`: An updated `NearestNeighbors` with new observation names.
#' }
#'
#' @export
#' @rdname NearestNeighbors-methods
setMethod(
  f = "names<-",
  signature = c("NearestNeighbors", "character"),
  definition = function(x, value) {
    if (length(value) == 0) {
      x@NAMES <- NULL
      return(x)
    }
    stopifnot(length(value) == nrow(x@nnidx))
    stopifnot(!anyDuplicated(value))
    x@NAMES <- value
    x
  }
)

#' @export
#' @rdname NearestNeighbors-methods
setMethod(
  f = "names<-",
  signature = c("NearestNeighbors", "NULL"),
  definition = function(x, value) {
    x@NAMES <- value
    x
  }
)

## dim #########################################################################

#' @returns
#' \itemize{
#' \item `dim`: A two-element vector represent the number of observations and
#' the number of nearest neighbors for each observation.
#' }
#'
#' @examples
#' dim(nn)
#'
#' @export
#' @rdname NearestNeighbors-methods
setMethod(
  f = "dim",
  signature = "NearestNeighbors",
  definition = function(x) dim(x@nnidx)
)

## dimnames ####################################################################

#' @returns
#' \itemize{
#' \item `dimnames`: The names of observations. Note that `colnames(x)` will
#' always be `NULL`.
#' }
#'
#' @examples
#' dimnames(nn)
#'
#' @export
#' @rdname NearestNeighbors-methods
setMethod(
  f = "dimnames",
  signature = "NearestNeighbors",
  definition = function(x) list(names(x), NULL)
)

#' @returns
#' \itemize{
#' \item `dimnames<-`: An updated `x` with new observation names. Note that
#' `colnames(x)<-` will be ignored.
#' }
#'
#' @examples
#' rownames(nn) <- sample(letters, nrow(nn))
#'
#' @export
#' @rdname NearestNeighbors-methods
setMethod(
  f = "dimnames<-",
  signature = c("NearestNeighbors", "list"),
  definition = function(x, value) {
    names(x) <- value[[1L]]
    x
  }
)

## parallel_slot_names #########################################################

#' @examples
#' nn[c(1, 3, 4, 5)]
#' nn[sample(names(nn), 10)]
#'
#' @export
#' @rdname NearestNeighbors-methods
setMethod(
  f = "parallel_slot_names",
  signature = "NearestNeighbors",
  definition = function(x) c("nnidx", "nndist", "NAMES", callNextMethod())
)

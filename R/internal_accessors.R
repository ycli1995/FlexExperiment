#' @importFrom SingleCellExperiment int_metadata
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomeInfoDb isCircular isCircular<- genome genome<- seqinfo
#' seqinfo<- seqlengths seqlengths<- seqlevels seqlevels<- seqnames seqnames<-
#' @importFrom S4Vectors DataFrame mcols mcols<- metadata metadata<-
#'
#' @importFrom easy.utils fastIntersect
#' @importFrom rlang is_character is_function
#' @importFrom BiocGenerics path updateObject
#' @importFrom methods as
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @importClassesFrom S4Vectors SimpleList
NULL

#' Helper functions for internal accessors
#'
#' Functions to access internal fields for any class that has implement methods
#' of `r .doc_links("int_colData")` and `r .doc_links("int_elementMetadata")`.
#'
#' @param x An object
#' @param ... `r .dot_param`
#'
#' @name internal-helpers
NULL

#' @param getfun A function to fetch the internal `r .doc_links("DataFrame")`.
#' Typically one of the following:
#' \itemize{
#' \item `r .doc_links("int_colData")`
#' \item `r .doc_links("int_elementMetadata")`
#' }
#' @param key Name of the fetched `r .doc_links("DataFrame")`. For example,
#' 'reducedDims' and 'colPairs'.
#'
#' @export
#' @rdname internal-helpers
getInternalAll <- function(x, getfun, key) {
  x <- updateObject(x)
  as(getfun(x)[[key]], "SimpleList")
}

#' @param value An object to set into the internal `r .doc_links("DataFrame")`.
#' @param setfun A function to set the updated `r .doc_links("DataFrame")` back
#' into `x`. #' Typically one of the following:
#' \itemize{
#' \item `r .doc_links("int_colData<-")`
#' \item `r .doc_links("int_elementMetadata<-")`
#' }
#' @param funstr A string specify the name of wrapper function that calls the
#' corresponding setter function.
#' @param name.prefix A string specify the prefix for unnamed elements. Default
#' is 'unnamed'.
#'
#' @importFrom S4Vectors I
#' @export
#' @rdname internal-helpers
setInternalAll <- function(
    x, value,
    getfun,
    setfun,
    key,
    prefix = NULL
) {
  x <- updateObject(x)
  int.data <- getfun(x)
  if (length(value) == 0L) {
    int.data[[key]] <- new2("DFrame", nrows = nrow(int.data), check = FALSE)
    return(setfun(x, int.data))
  }
  names(value) <- .clean_internal_names(
    names(value),
    N = length(value),
    what = "names(value)",
    prefix = prefix
  )
  collected <- do.call(
    what = DataFrame,
    args = c(lapply(value, I), list(row.names = NULL, check.names = FALSE))
  )
  if (is(value, "Annotated")) {
    metadata(collected) <- metadata(value)
  }
  if (is(value, "Vector")) {
    mcols(collected) <- mcols(value)
  }
  int.data[[key]] <- collected
  setfun(x, int.data)
}

#' @export
#' @rdname internal-helpers
getInternalNames <- function(x, getfun, key) {
  x <- updateObject(x)
  colnames(getfun(x)[[key]])
}

#' @export
#' @rdname internal-helpers
setInternalNames <- function(
    x,
    value,
    getfun,
    setfun,
    key,
    funstr,
    prefix = NULL
) {
  x <- updateObject(x)
  int.data <- getfun(x)
  value <- .clean_internal_names(
    value,
    N = ncol(int.data[[key]]),
    what = sprintf("%s(x)", funstr),
    prefix = prefix
  )
  colnames(int.data[[key]]) <- value
  setfun(x, int.data)
}

#' @param index Which element to fetch or set the internal data.
#' @param substr A string specify the parameter name of `index`. For example,
#' 'type' for `r .doc_links("reducedDim")` or `r .doc_links("colPair")`
#' @param type A string specify the type of `index`. For example, 'character' or
#' 'numeric'.
#'
#' @export
#' @rdname internal-helpers
getInternalOne <- function(x, i, getfun, key, funstr) {
  x <- updateObject(x)
  int.data <- getfun(x)[[key]]
  tryCatch(int.data[, i], error = function(err) {
    e <- sprintf("invalid subscript '%s' in '%s': \n", as.character(i), funstr)
    stop(e, conditionMessage(err))
  })
}

#' @export
#' @rdname internal-helpers
setInternalOne <- function(
    x, i,
    value,
    getfun,
    setfun,
    key,
    namestr,
    substr,
    funstr,
    prefix = NULL
) {
  i <- i[1]
  x <- updateObject(x)
  int.data <- getfun(x)
  if (!is.character(i)) {
    if (i > ncol(int.data[[key]])) {
      fmt <- "'%s' out of bounds in '%s': %s"
      stop(sprintf(fmt, substr, funstr, as.character(i)))
    }
  }
  int.data[[key]][[i]] <- value
  names(int.data[[key]]) <- .clean_internal_names(
    names(int.data[[key]]),
    N = ncol(int.data[[key]]),
    what = sprintf("%s(x)", namestr),
    prefix = prefix
  )
  setfun(x, int.data)
}

getInternalMissing <- function(x, basefun, namefun, substr, ...) {
  if (length(namefun(x)) == 0) {
    funstr <- as.character(substitute(basefun))
    fmt <- "'%s(<%s>, %s=\"missing\")' failed: no avaliable %s"
    stop(sprintf(fmt, funstr, class(x)[1], substr, funstr))
  }
  basefun(x, 1L, ...)
}

setInternalMissing <- function(x, value, basefun, namefun, ...) {
  if (length(namefun(x)) > 0) {
    i <- 1L
  } else {
    i <- paste0(.unnamed, 1L)
  }
  basefun(x, i, ..., value = value)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.unnamed <- "unnamed"

.clean_internal_names <- function(names, N, what, prefix = NULL) {
  prefix <- prefix %||% "unnamed"
  if (length(names) == 0 && N > 0) {
    fmt <- "%s is empty, replacing with '%s'"
    warning(sprintf(fmt, what, prefix), immediate. = TRUE, call. = FALSE)
    return(paste0(prefix, seq_len(N)))
  }
  if (!any(empty <- (names == "") | is.na(names))) {
    return(names)
  }
  fmt <- "%s contains empty strings, replacing with prefix '%s'"
  warning(sprintf(fmt, what, prefix), immediate. = TRUE, call. = FALSE)
  names[empty] <- paste0(prefix, which(empty))
  names
}

.check_equal_for_setter <- function(x, y, xname, yname, funstr) {
  fmt <- "%s failed:\n '%s' (%d) should be the same as '%s' (%d)"
  err <- sprintf(fmt, funstr, xname, x, yname, y)
  checkEqual(x, y, err = err)
}

.check_identical_for_setter <- function(x, y, xname, yname, funstr) {
  fmt <- "%s failed:\n '%s' should be the same as '%s'"
  err <- sprintf(fmt, funstr, xname, yname)
  checkIdentical(x, y, err = err)
}



#' @include utils.R
#' @include FlexExperiment-class.R
#' @include internal_accessors.R
#'
#' @importFrom SingleCellExperiment counts counts<- cpm cpm<- logcounts
#' logcounts<- normcounts normcounts<- tpm tpm<-
#' @importFrom IRanges relist PartitioningByEnd
#' @importFrom S4Vectors make_zero_col_DFrame setValidity2
#' @importFrom methods setClass
NULL

#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames
#' reducedDimNames<- reducedDims reducedDims<-
.red_key <- "reducedDims"

#' @importFrom SingleCellExperiment colPair colPair<- colPairNames
#' colPairNames<- colPairs colPairs<-
.colp_key <- "colPairs"

#' @importFrom SingleCellExperiment rowPair rowPair<- rowPairNames
#' rowPairNames<- rowPairs rowPairs<-
.rowp_key <- "rowPairs"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class ########################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The SingleCellFlexExperiment
#'
#' The `SingleCellFlexExperiment` class is designed to represent single-cell
#' data. It inherits from \code{\link{LayeredExperiment}} class and is used in
#' the same manner. By default, this class supplies storage of dimensionality
#' reduction results (\code{\link{reducedDims}}), row or column pairing data
#' (\code{\link{rowPairs}} and \code{\link{colPairs}}) and column nearest
#' neighbor (NN) results (\code{\link{colNeighbors}}).
#'
#' @name SingleCellFlexExperiment
#' @docType class
#' @exportClass SingleCellFlexExperiment
setClass("SingleCellFlexExperiment", contains = "FlexExperiment", slots = c(
  int_colData = "DataFrame",
  int_elementMetadata = "DataFrame",
  int_metadata = "list"
))

setMethod("initialize", "SingleCellFlexExperiment", function(.Object, ...) {
  .Object@int_metadata <- c(
    .Object@int_metadata,
    list(version = packageVersion(.packageName))
  )
  callNextMethod(.Object, ...)
})

## Constructor #################################################################

#' @param assays A \code{\link{list}} or `r .doc_links("SimpleList")` of
#' matrix-like objects. Also supports a pre-built \code{\link{LayeredAssays}}.
#' Names of rows and columns must present in those matrix-like objects.
#' @param rowData A `r .doc_links("DataFrame")` object summarizing the rows of
#' `assays`. Row names must present to match each row to corresponding rows in
#' `assays`. If `rowData` is not provided, the output row names will use the
#' union of row names in `assays`.
#' @param rowRanges A `r .doc_links("GRanges")` or `r .doc_links("GRangesList")`
#' object describing the ranges of interest. The length of `rowRanges` must be
#' equal to the row number of `assays`.
#' @param colData A `r .doc_links("DataFrame")` object summarizing the columns of
#' `assays`. Row names must present to match each row to corresponding columns
#' in `assays`. If `colData` is not provided, the output column names will use
#' the union of column names in `assays`.
#' @param reducedDims A list of any number of matrix-like objects representing
#' dimensionality reduction results, each of which should have the same number
#' of rows as the output `SingleCellFlexExperiment`.
#' @param rowPairs A list of any number of `r .doc_links("CsparseMatrix")`
#' representing the pairing data of rows, each of which should inherit from
#' `dMatrix` and contain the same numbers of rows and columns, matching rows of
#' the output `SingleCellFlexExperiment`.
#' @param colPairs A list of any number of `r .doc_links("CsparseMatrix")`
#' representing the pairing data of columns, each of which should inherit from
#' `dMatrix` and contain the same numbers of rows and columns, matching columns
#' of the output `SingleCellFlexExperiment`.
#' @param colNeighbors A list of any number of \code{\link{NearestNeighbors}}
#' objects, each of which should have the same number of observations as the
#' column number of the output `SingleCellFlexExperiment`.
#' @param metadata An optional \code{\link{list}} of arbitrary content
#' describing the overall experiment.
#'
#' @importFrom GenomicRanges GRangesList
#'
#' @export
#' @rdname SingleCellFlexExperiment
SingleCellFlexExperiment <- function(
    ...,
    reducedDims = list(),
    rowPairs = list(),
    colPairs = list()
) {
  se <- FlexExperiment(...)
  .FE_to_SCFE(
    se,
    reducedDims = reducedDims,
    rowPairs = rowPairs,
    colPairs = colPairs
  )
}

.FE_to_SCFE <- function(
    fe,
    reducedDims = list(),
    rowPairs = list(),
    colPairs = list()
) {
  int_elementMetadata <- .empty_int_DFrame(nrow(fe), .rowp_key)
  int_colData <- .empty_int_DFrame(ncol(fe), c(.red_key, .colp_key))
  out <- new(
    "SingleCellFlexExperiment", fe,
    int_colData = int_colData,
    int_elementMetadata = int_elementMetadata
  )
  reducedDims(out) <- reducedDims
  rowPairs(out) <- rowPairs
  colPairs(out) <- colPairs
  return(out)
}

.empty_int_DFrame <- function(nrow, keys = NULL) {
  df <- new2("DFrame", nrows = nrow, check = FALSE)
  for (i in keys) {
    df[[i]] <- new2("DFrame", nrows = nrow, check = FALSE)
  }
  df
}

## Validity ####################################################################

validSCFEIntData <- function(x, getfun, dimfun, keys, immediate. = TRUE) {
  intdata <- getfun(x)
  funstr <- as.character(substitute(getfun))
  err <- NULL
  if (nrow(intdata) != dimfun(x)) {
    dimstr <- as.character(substitute(dimfun))
    e <- sprintf("'nrow' of '%s' not equal to '%s(x)'", funstr, dimstr)
    err <- getErrors(e, err, immediate. = immediate.)
  }
  if (any(!keys %in% names(intdata))) {
    miss.keys <- setdiff(keys, names(intdata))
    fmt <- "Required key(s) not found in '%s(x)': %s"
    e <- sprintf(fmt, funstr, paste(miss.keys, collapse = ", "))
    err <- getErrors(e, err, immediate.)
  }
  invisible(err)
}

setValidity2("SingleCellFlexExperiment", function(x) {
  invisible(c(
    validFE(x, immediate. = FALSE),
    validSCFEIntData(x, int_colData, ncol, c(.red_key, .colp_key), FALSE),
    validSCFEIntData(x, int_elementMetadata, nrow, .rowp_key, FALSE)
  ))
})

#' @export
setMethod(
  f = "updateObject",
  signature = "SingleCellFlexExperiment",
  definition = function(object, ..., verbose = FALSE) {
    object@int_colData <- updateObject(
      object@int_colData, ...,
      verbose = verbose
    )
    object@int_elementMetadata <- updateObject(
      object@int_elementMetadata, ...,
      verbose = verbose
    )
    int_metadata(object)$version <- NULL
    object@int_metadata <- updateObject(
      object@int_metadata, ...,
      verbose = verbose
    )
    int_metadata(object)$version <- packageVersion(.packageName)
    callNextMethod()
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SingleCellFlexExperiment-methods #############################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Methods for SingleCellFlexExperiment
#'
#' Accessor methods specific for \code{\link{SingleCellFlexExperiment}}.
#'
#' @param x,object A \code{\link{SingleCellFlexExperiment}} object.
#' @param ... `r .dot_param`
#'
#' @seealso
#' S4 generics defined in \pkg{SingleCellExperiment}:
#' \itemize{
#' \item `r .doc_links("reducedDims")`
#' \item `r .doc_links("colPairs")`
#' \item `r .doc_links("rowPairs")`
#' }
#' S4 generics defined in \pkg{OmicsArtifacts}:
#' \itemize{
#' \item \code{\link{colNeighbors}}
#' }
#'
#' @name SingleCellFlexExperiment-methods
NULL

## show ########################################################################

#' @export
showSCFE <- function(object) {
  showFE(object)
  coolcat("reducedDimNames(%d): %s\n", reducedDimNames(object))
  coolcat("colPairs(%d): %s\n", colPairNames(object))
  coolcat("rowPairs(%d): %s\n", rowPairNames(object))
}

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod("show", "SingleCellFlexExperiment", showSCFE)

setMethod("int_colData", "SingleCellFlexExperiment", function(x) x@int_colData)

setMethod("int_colData<-", "SingleCellFlexExperiment", function(x, value) {
  x@int_colData <- value
  validSCFEIntData(x, int_colData, ncol, c(.red_key, .colp_key))
  x
})

setMethod("int_elementMetadata", "SingleCellFlexExperiment", function(x) {
  x@int_elementMetadata
})

setMethod(
  f = "int_elementMetadata<-",
  signature = "SingleCellFlexExperiment",
  definition = function(x, value) {
    x@int_elementMetadata <- value
    validSCFEIntData(x, int_elementMetadata, nrow, .rowp_key)
    x
  }
)

setMethod("int_metadata", "SingleCellFlexExperiment", function(x) {
  x@int_metadata
})

setMethod("int_metadata<-", "SingleCellFlexExperiment", function(x, value) {
  x@int_metadata <- value
  x
})

## reducedDims #################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDims",
  signature = "SingleCellFlexExperiment",
  definition = function(x, withDimnames = TRUE, ...) {
    value <- getInternalAll(x, getfun = int_colData, key = .red_key)
    if (!withDimnames) {
      return(value)
    }
    for (i in seq_along(value)) {
      rownames(value[[i]]) <- colnames(x)
    }
    value
  }
)

## reducedDims<- ###############################################################

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

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDims<-",
  signature = "SingleCellFlexExperiment",
  definition = function(x, withDimnames = TRUE, ..., value) {
    funstr <- sprintf("`reducedDims<-`(<%s>, ...):\n ", class(x)[1])
    if (withDimnames) {
      for (i in seq_along(value)) {
        .check_identical_for_setter(
          x = rownames(value[[i]]),
          y = colnames(x),
          xname = paste0("rownames(value[[", i, "]])"),
          yname = "colnames(x)",
          funstr = funstr
        )
      }
    } else {
      for (i in seq_along(value)) {
        .check_equal_for_setter(
          x = nrow(value[[i]]),
          y = ncol(x),
          xname = paste0("nrow(value[[", i, "]])"),
          yname = "ncol(x)",
          funstr = funstr
        )
      }
    }
    setInternalAll(
      x = x,
      value = value,
      getfun = int_colData,
      setfun = `int_colData<-`,
      key = .red_key,
      ...
    )
  }
)

## reducedDimNames #############################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod("reducedDimNames", "SingleCellFlexExperiment", function(x) {
  getInternalNames(x, getfun = int_colData, key = .red_key)
})

## reducedDimNames<- ###########################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDimNames<-",
  signature = c("SingleCellFlexExperiment", "character"),
  definition = function(x, value) {
    setInternalNames(
      x, value,
      getfun = int_colData,
      setfun = `int_colData<-`,
      key = .red_key,
      funstr = "reducedDimNames"
    )
  }
)

## reducedDim ##################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDim",
  signature = c("SingleCellFlexExperiment", "missing"),
  definition = function(x, type, withDimnames = TRUE, ...) {
    getInternalMissing(
      x,
      basefun = reducedDim,
      namefun = reducedDimNames,
      substr = "type",
      withDimnames = withDimnames,
      ...
    )
  }
)

.get_reducedDim <- function(x, i, type, withDimnames = TRUE, ...) {
  funstr <- sprintf("reducedDim(<%s>, type = <%s>", class(x)[1], type)
  out <- getInternalOne(
    x, i,
    getfun = int_colData,
    key = .red_key,
    funstr = funstr
  )
  if (withDimnames) {
    rownames(out) <- colnames(x)
  }
  out
}

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDim",
  signature = c("SingleCellFlexExperiment", "integer"),
  definition = function(x, type, withDimnames = TRUE, ...) {
    .get_reducedDim(
      x = x,
      i = type,
      type = "integer",
      withDimnames = withDimnames,
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDim",
  signature = c("SingleCellFlexExperiment", "numeric"),
  definition = function(x, type, withDimnames = TRUE, ...) {
    .get_reducedDim(
      x = x,
      i = as.integer(type),
      type = "numeric",
      withDimnames = withDimnames,
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDim",
  signature = c("SingleCellFlexExperiment", "character"),
  definition = function(x, type, withDimnames = TRUE, ...) {
    .get_reducedDim(
      x = x,
      i = type,
      type = "character",
      withDimnames = withDimnames,
      ...
    )
  }
)

## reducedDim<- ################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDim<-",
  signature = c("SingleCellFlexExperiment", "missing"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    setInternalMissing(
      x, value,
      basefun = `reducedDim<-`,
      namefun = reducedDimNames,
      withDimnames = withDimnames,
      ...
    )
  }
)

.set_reducedDim <- function(x, i, type, value, withDimnames = TRUE, ...) {
  funstr <- sprintf('`reducedDim<-`(<%s>, type = <%s>, ...)', class(x)[1], type)
  if (withDimnames) {
    .check_identical_for_setter(
      x = rownames(value),
      y = colnames(x),
      xname = "rownames(value)",
      yname = "colnames(x)",
      funstr = funstr
    )
  } else {
    .check_equal_for_setter(
      x = nrow(value),
      y = ncol(x),
      xname = "nrow(value)",
      yname = "ncol(x)",
      funstr = funstr
    )
  }
  setInternalOne(
    x, i, value,
    getfun = int_colData,
    setfun = `int_colData<-`,
    key = .red_key,
    funstr = funstr,
    namestr = "reducedDimNames",
    substr = "type"
  )
}

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDim<-",
  signature = c("SingleCellFlexExperiment", "integer"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    .set_reducedDim(
      x = x,
      i = type,
      type = "integer",
      value = value,
      withDimnames = withDimnames,
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDim<-",
  signature = c("SingleCellFlexExperiment", "numeric"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    .set_reducedDim(
      x = x,
      i = as.integer(type),
      type = "numeric",
      value = value,
      withDimnames = withDimnames,
      ...
    )
  }
)

#' @importFrom SingleCellExperiment reducedDim<-
#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "reducedDim<-",
  signature = c("SingleCellFlexExperiment", "character"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    .set_reducedDim(
      x = x,
      i = type,
      type = "character",
      value = value,
      withDimnames = withDimnames,
      ...
    )
  }
)

## colPairs ####################################################################

.get_dimPairs <- function(x, getfun, dimfun, key, asSparse = FALSE, ...) {
  value <- getInternalAll(x, getfun = getfun, key = key)
  if (asSparse) {
    for (i in seq_along(value)) {
      names(value[[i]]) <- dimfun(x)
      value[[i]] <- hitsToMat(value[[i]])
    }
    return(value)
  }
  for (i in seq_along(value)) {
    value[[i]] <- value[[i]]@hits
  }
  value
}

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPairs",
  signature = "SingleCellFlexExperiment",
  definition = function(x, asSparse = FALSE, ...) {
    .get_dimPairs(x, int_colData, colnames, .colp_key, asSparse, ...)
  }
)

## colPairs<- ##################################################################

.set_dimPairs <- function(
    x, value,
    getfun, setfun,
    dimfun, nfun,
    funstr, key,
    ...
) {
  for (i in seq_along(value)) {
    if (!inherits(value[[i]], "PairedHits")) {
      value[[i]] <- as(value[[i]], "PairedHits")
    }
    if (length(names(value[[i]])) > 0) {
      .check_identical_for_setter(
        x = names(value[[i]]),
        y = dimfun(x),
        xname = paste0("names(value[[", i, "]])"),
        yname = sprintf("%s(x)", as.character(substitute(dimfun))),
        funstr = funstr
      )
    } else {
      .check_equal_for_setter(
        x = length(value[[i]]),
        y = nfun(x),
        xname = paste0("length(value[[", i, "]])"),
        yname = sprintf("%s(x)", as.character(substitute(nfun))),
        funstr = funstr
      )
    }
  }
  setInternalAll(x, value, getfun, setfun, key = key, ...)
}

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod("colPairs<-", "SingleCellFlexExperiment", function(x, value) {
  funstr <- sprintf("`colPairs<-`(<%s>, ...):\n ", class(x)[1])
  .set_dimPairs(
    x = x,
    value = value,
    getfun = int_colData,
    setfun = `int_colData<-`,
    dimfun = colnames,
    nfun = ncol,
    funstr = funstr,
    key = .colp_key
  )
})

## colPairNames ################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod("colPairNames", "SingleCellFlexExperiment", function(x) {
  getInternalNames(x, getfun = int_colData, key = .colp_key)
})

## colPairNames<- ##############################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPairNames<-",
  signature = c("SingleCellFlexExperiment", "character"),
  definition = function(x, value) {
    setInternalNames(
      x, value,
      getfun = int_colData,
      setfun = `int_colData<-`,
      key = .colp_key,
      funstr = "colPairNames"
    )
  }
)

## colPair #####################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPair",
  signature = c("SingleCellFlexExperiment", "missing"),
  definition = function(x, type, asSparse = FALSE, ...) {
    getInternalMissing(
      x = x,
      basefun = colPair,
      namefun = colPairNames,
      substr = "type",
      asSparse = asSparse
    )
  }
)

.get_dimPair <- function(
    x, i,
    getfun, dimfun,
    key, funstr,
    asSparse = FALSE,
    ...
) {
  out <- getInternalOne(
    x, i,
    getfun = getfun,
    key = key,
    funstr = funstr
  )
  if (asSparse) {
    names(out) <- dimfun(x)
    return(hitsToMat(out, ...))
  }
  out@hits
}

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPair",
  signature = c("SingleCellFlexExperiment", "integer"),
  definition = function(x, type, asSparse = FALSE, ...) {
    funstr <- sprintf("colPair(<%s>, type = <%s>)", class(x)[1], "integer")
    .get_dimPair(
      x, type,
      getfun = int_colData,
      dimfun = colnames,
      key = .colp_key,
      funstr = funstr,
      asSparse = asSparse,
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPair",
  signature = c("SingleCellFlexExperiment", "numeric"),
  definition = function(x, type, asSparse = FALSE, ...) {
    funstr <- sprintf("colPair(<%s>, type = <%s>)", class(x)[1], "numeric")
    .get_dimPair(
      x, as.integer(type),
      getfun = int_colData,
      dimfun = colnames,
      key = .colp_key,
      funstr = funstr,
      asSparse = asSparse,
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPair",
  signature = c("SingleCellFlexExperiment", "character"),
  definition = function(x, type, asSparse = FALSE, ...) {
    funstr <- sprintf("colPair(<%s>, type = <%s>)", class(x)[1], "character")
    .get_dimPair(
      x, type,
      getfun = int_colData,
      dimfun = colnames,
      key = .colp_key,
      funstr = funstr,
      asSparse = asSparse,
      ...
    )
  }
)

## colPair<- ###################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPair<-",
  signature = c("SingleCellFlexExperiment", "missing"),
  definition = function(x, type, ..., value) {
    setInternalMissing(
      x, value,
      basefun = `colPair<-`,
      namefun = colPairNames,
      ...
    )
  }
)

.set_dimPair <- function(
    x, i, value,
    getfun, setfun,
    dimfun, nfun,
    key, funstr, namestr,
    ...
) {
  if (!inherits(value, "PairedHits")) {
    value <- as(value, "PairedHits")
  }
  if (length(names(value)) > 0) {
    .check_identical_for_setter(
      x = names(value),
      y = dimfun(x),
      xname = "names(value)",
      yname = sprintf("%s(x)", as.character(substitute(dimfun))),
      funstr = funstr
    )
  } else {
    .check_equal_for_setter(
      x = length(value),
      y = nfun(x),
      xname = "length(value)",
      yname = sprintf("%s(x)", as.character(substitute(nfun))),
      funstr = funstr
    )
  }
  setInternalOne(
    x, i, value,
    getfun = getfun,
    setfun = setfun,
    key = key,
    funstr = funstr,
    namestr = namestr,
    substr = "type",
    ...
  )
}

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPair<-",
  signature = c("SingleCellFlexExperiment", "integer"),
  definition = function(x, type, ..., value) {
    fstr <- sprintf("`colPair<-`(<%s>, type = <%s>)", class(x)[1], "integer")
    .set_dimPair(
      x, type, value,
      getfun = int_colData,
      setfun = `int_colData<-`,
      dimfun = colnames,
      nfun = ncol,
      key = .colp_key,
      funstr = fstr,
      namestr = "colPairNames",
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPair<-",
  signature = c("SingleCellFlexExperiment", "numeric"),
  definition = function(x, type, ..., value) {
    fstr <- sprintf("`colPair<-`(<%s>, type = <%s>)", class(x)[1], "numeric")
    .set_dimPair(
      x, as.integer(type), value,
      getfun = int_colData,
      setfun = `int_colData<-`,
      dimfun = colnames,
      nfun = ncol,
      key = .colp_key,
      funstr = fstr,
      namestr = "colPairNames",
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "colPair<-",
  signature = c("SingleCellFlexExperiment", "character"),
  definition = function(x, type, ..., value) {
    fstr <- sprintf("`colPair<-`(<%s>, type = <%s>)", class(x)[1], "character")
    .set_dimPair(
      x, type, value,
      getfun = int_colData,
      setfun = `int_colData<-`,
      dimfun = colnames,
      nfun = ncol,
      key = .colp_key,
      funstr = fstr,
      namestr = "colPairNames",
      ...
    )
  }
)

## rowPairs ####################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPairs",
  signature = "SingleCellFlexExperiment",
  definition = function(x, asSparse = FALSE, ...) {
    .get_dimPairs(x, int_elementMetadata, rownames, .rowp_key, asSparse, ...)
  }
)

## rowPairs<- ##################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod("rowPairs<-", "SingleCellFlexExperiment", function(x, value) {
  funstr <- sprintf("`rowPairs<-`(<%s>, ...):\n ", class(x)[1])
  .set_dimPairs(
    x = x,
    value = value,
    getfun = int_elementMetadata,
    setfun = `int_elementMetadata<-`,
    dimfun = rownames,
    nfun = nrow,
    funstr = funstr,
    key = .rowp_key
  )
})

## rowPairNames ################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod("rowPairNames", "SingleCellFlexExperiment", function(x) {
  getInternalNames(x, getfun = int_elementMetadata, key = .rowp_key)
})

## rowPairNames<- ##############################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPairNames<-",
  signature = c("SingleCellFlexExperiment", "character"),
  definition = function(x, value) {
    setInternalNames(
      x, value,
      getfun = int_elementMetadata,
      setfun = `int_elementMetadata<-`,
      key = .rowp_key,
      funstr = "rowPairNames"
    )
  }
)

## rowPair #####################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPair",
  signature = c("SingleCellFlexExperiment", "missing"),
  definition = function(x, type, asSparse = FALSE, ...) {
    getInternalMissing(
      x = x,
      basefun = rowPair,
      namefun = rowPairNames,
      substr = "type",
      asSparse = asSparse
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPair",
  signature = c("SingleCellFlexExperiment", "integer"),
  definition = function(x, type, asSparse = FALSE, ...) {
    funstr <- sprintf("rowPair(<%s>, type = <%s>)", class(x)[1], "integer")
    .get_dimPair(
      x, type,
      getfun = int_elementMetadata,
      dimfun = rownames,
      key = .rowp_key,
      funstr = funstr,
      asSparse = asSparse,
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPair",
  signature = c("SingleCellFlexExperiment", "numeric"),
  definition = function(x, type, asSparse = FALSE, ...) {
    funstr <- sprintf("rowPair(<%s>, type = <%s>)", class(x)[1], "numeric")
    .get_dimPair(
      x, as.integer(type),
      getfun = int_elementMetadata,
      dimfun = rownames,
      key = .rowp_key,
      funstr = funstr,
      asSparse = asSparse,
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPair",
  signature = c("SingleCellFlexExperiment", "character"),
  definition = function(x, type, asSparse = FALSE, ...) {
    funstr <- sprintf("rowPair(<%s>, type = <%s>)", class(x)[1], "character")
    .get_dimPair(
      x, type,
      getfun = int_elementMetadata,
      dimfun = rownames,
      key = .rowp_key,
      funstr = funstr,
      asSparse = asSparse,
      ...
    )
  }
)

## rowPair<- ###################################################################

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPair<-",
  signature = c("SingleCellFlexExperiment", "missing"),
  definition = function(x, type, ..., value) {
    setInternalMissing(
      x, value,
      basefun = `rowPair<-`,
      namefun = rowPairNames,
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPair<-",
  signature = c("SingleCellFlexExperiment", "integer"),
  definition = function(x, type, ..., value) {
    fstr <- sprintf("`rowPair<-`(<%s>, type = <%s>)", class(x)[1], "integer")
    .set_dimPair(
      x, type, value,
      getfun = int_elementMetadata,
      setfun = `int_elementMetadata<-`,
      dimfun = rownames,
      nfun = nrow,
      key = .rowp_key,
      funstr = fstr,
      namestr = "rowPairNames",
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPair<-",
  signature = c("SingleCellFlexExperiment", "numeric"),
  definition = function(x, type, ..., value) {
    fstr <- sprintf("`rowPair<-`(<%s>, type = <%s>)", class(x)[1], "numeric")
    .set_dimPair(
      x, as.integer(type), value,
      getfun = int_elementMetadata,
      setfun = `int_elementMetadata<-`,
      dimfun = rownames,
      nfun = nrow,
      key = .rowp_key,
      funstr = fstr,
      namestr = "rowPairNames",
      ...
    )
  }
)

#' @export
#' @rdname SingleCellFlexExperiment-methods
setMethod(
  f = "rowPair<-",
  signature = c("SingleCellFlexExperiment", "character"),
  definition = function(x, type, ..., value) {
    fstr <- sprintf("`rowPair<-`(<%s>, type = <%s>)", class(x)[1], "character")
    .set_dimPair(
      x, type, value,
      getfun = int_elementMetadata,
      setfun = `int_elementMetadata<-`,
      dimfun = rownames,
      nfun = nrow,
      key = .rowp_key,
      funstr = fstr,
      namestr = "rowPairNames",
      ...
    )
  }
)

#' @export
#' @rdname FlexExperiment-subset
setMethod(
  f = "[",
  signature = "SingleCellFlexExperiment",
  definition = function(x, i, j, ..., drop = FALSE) {
    x <- updateObject(x)
    if (!missing(i)) {
      ii <- .subset_names2idx(i, rownames(x))
      x@int_elementMetadata <- int_elementMetadata(x)[ii, , drop = FALSE]
    }
    if (!missing(j)) {
      jj <- .subset_names2idx(j, colnames(x))
      x@int_colData <- int_colData(x)[jj, , drop = FALSE]
    }
    callNextMethod()
  }
)

GET_ASSAY <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ...) assay(object, i = exprs_values, ...)
}

SET_ASSAY <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ..., value) {
    assay(object, i = exprs_values, ...) <- value
    object
  }
}

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("counts", "SingleCellFlexExperiment", GET_ASSAY("counts"))

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("counts<-", "SingleCellFlexExperiment", SET_ASSAY("counts"))

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("logcounts", "SingleCellFlexExperiment", GET_ASSAY("logcounts"))

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("logcounts<-", "SingleCellFlexExperiment", SET_ASSAY("logcounts"))

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("normcounts", "SingleCellFlexExperiment", GET_ASSAY("normcounts"))

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("normcounts<-", "SingleCellFlexExperiment", SET_ASSAY("normcounts"))

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("cpm", "SingleCellFlexExperiment", GET_ASSAY("cpm"))

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("cpm<-", "SingleCellFlexExperiment", SET_ASSAY("cpm"))

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("tpm", "SingleCellFlexExperiment", GET_ASSAY("tpm"))

#' @export
#' @rdname SingleCellFlexExperiment-assays
setMethod("tpm<-", "SingleCellFlexExperiment", SET_ASSAY("tpm"))




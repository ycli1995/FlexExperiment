#' @include utils.R
#'
#' @importFrom FlexAssays assayClasses autoDimnames colFlex FlexAssays getAssays
#' getErrors getOneAssay rowFlex setOneAssay showS4Title
#' @importFrom SummarizedExperiment assay assay<- assayNames assayNames<- assays
#' assays<- colData colData<- rowData rowData<- rowRanges rowRanges<-
#' @importFrom SingleCellExperiment int_colData int_colData<-
#' int_elementMetadata int_elementMetadata<- int_metadata int_metadata<-
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom S4Vectors coolcat DataFrame extractROWS mcols mcols<- metadata
#' metadata<- new2 setValidity2
#' @importFrom IRanges PartitioningByEnd relist
#' @importFrom BiocGenerics colnames<- rownames<- updateObject
#' @importFrom methods slot<-
#' @importFrom utils .DollarNames
#'
#' @importClassesFrom GenomicRanges GenomicRanges_OR_GRangesList
#' @importClassesFrom S4Vectors Annotated DataFrame DFrame
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class ########################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The FlexExperiment class
#'
#' The `FlexExperiment` class is a container where matrix-like objects are
#' stored in a \code{\link{FlexAssays}}. The rows represent features of
#' interest (e.g. genes, etc...) and columns represent samples (cells).
#'
#' @slot assays A \code{\link{FlexAssays}} object containing the layered
#' matrix-like objects.
#' @slot rowData,colData `r .doc_links("DataFrame")` that summarize the metadata
#' of rows/columns. The row number of `row(col)Data` must match the rows/columns
#' in `assays`.
#' @slot int_colData,int_elementMetadata `r .doc_links("DataFrame")` of internal
#' column/row metadata, with number of rows equal to `dim(x)`.
#' @slot int_metadata A list of internal metadata
#'
#' @details
#' In `FlexExperiment`, both row names and column names are only stored in
#' `assays`.
#'
#' @name FlexExperiment
#' @docType class
#'
#' @exportClass FlexExperiment
#' @rdname FlexExperiment
setClass("FlexExperiment", contains = "Annotated", slots = c(
  rowRanges = "GenomicRanges_OR_GRangesList",
  colData = "DataFrame",
  assays = "FlexAssays"
))

## Constructor #################################################################

#' @export
FlexExperiment <- function(
    assays = List(),
    rowData = NULL,
    rowRanges = GRangesList(),
    colData = DataFrame(),
    metadata = list(),
    assayClasses = NULL,
    rowFlex = c("bounded", "fixed"),
    colFlex = c("bounded", "fixed")
) {
  rowFlex <- match.arg(rowFlex)
  colFlex <- match.arg(colFlex)
  assayClasses <- as.character(assayClasses)

  rnames <- NULL
  cnames <- NULL
  if (!is(colData, "DataFrame")) {
    colData <- as(colData, "DataFrame")
  }
  if (nrow(colData) > 0) {
    stopifnot(length(rownames(colData)) > 0)
    cnames <- rownames(colData)
  }

  if (length(rowRanges) > 0) {
    if (!is.null(rowData)) {
      stop("only one of 'rowData' and 'rowRanges' can be specified")
    }
    rowData <- mcols(rowRanges)
  }
  if (!is.null(rowData)) {
    if (is(rowData, "GenomicRanges_OR_GRangesList")) {
      rowRanges <- rowData
      rowData <- mcols(rowRanges)
    }
    if (!is(rowData, "DataFrame")) {
      rowData <- as(rowData, "DataFrame")
    }
  } else {
    rowData <- new2("DFrame", nrows = 0L, check = FALSE)
  }
  if (nrow(rowData) > 0) {
    stopifnot(length(rownames(rowData)) > 0)
    rnames <- rownames(rowData)
  }

  if (length(assays) > 0) {
    assays <- FlexAssays(
      assays,
      assayClasses = assayClasses,
      rowFlex = rowFlex,
      colFlex = colFlex
    )
    if (nrow(rowData) > 0) {
      .check_FE_dim(assays, rowData, rownames)
    }
    if (ncol(colData) > 0) {
      .check_FE_dim(assays, colData, colnames)
    }
  } else {
    assays <- FlexAssays(
      rownames = rnames,
      colnames = cnames,
      assayClasses = assayClasses,
      rowFlex = rowFlex,
      colFlex = colFlex
    )
  }
  if (nrow(rowData) == 0) {
    rowData <- new2("DFrame", nrows = nrow(assays), check = FALSE)
  }
  if (length(rowRanges) == 0) {
    rowRanges <- .init_rowRanges(rowData)
  }
  if (length(rowRanges) > 0) {
    names(rowRanges) <- NULL
  }
  if (nrow(colData) == 0) {
    colData <- new2("DFrame", nrows = ncol(assays), check = FALSE)
  } else {
    rownames(colData) <- NULL
  }
  new2(
    Class = "FlexExperiment",
    assays = assays,
    rowRanges = rowRanges,
    colData = colData,
    metadata = metadata,
    check = FALSE
  )
}

.init_rowRanges <- function(rowData) {
  partitioning <- PartitioningByEnd(
    x = integer(nrow(rowData)),
    names = rownames(rowData)
  )
  rowRanges <- relist(GRanges(), partitioning)
  mcols(rowRanges) <- rowData
  rowRanges
}

.check_FE_dim <- function(assays, dimData, dimfun) {
  if (.identical_fmatch(dimfun(assays), rownames(dimData))) {
    return(invisible(NULL))
  }
  funstr <- as.character(substitute(dimfun))
  datastr <- sub("names", "Data", funstr)
  stop(sprintf("%s(assays) differs from rownames(%s).", funstr, datastr))
}

## Validity ####################################################################

.valid_FE_rowRanges <- function(assays, rowRanges, immediate. = TRUE) {
  if (length(rowRanges) == nrow(assays)) {
    return(invisible(NULL))
  }
  fmt <- "length of 'rowRanges' (%d) must be equal nrow(assays) (%d)"
  err <- sprintf(fmt, length(rowRanges), nrow(assays))
  getErrors(err, immediate. = immediate.)
}

.valid_FE_ndim <- function(assays, dimData, dimfun, immediate. = TRUE) {
  if (dimfun(assays) == nrow(dimData)) {
    return(invisible(NULL))
  }
  funstr <- as.character(substitute(dimfun))
  datastr <- sub("n", "", funstr)
  fmt <- "%s(assays) (%d) differs from nrow(%sData) (%d)."
  getErrors(
    sprintf(fmt, funstr, dimfun(assays), datastr, nrow(dimData)),
    immediate. = immediate.
  )
}

validFE <- function(x, immediate. = TRUE) {
  invisible(c(
    .valid_FE_ndim(x@assays, x@colData, ncol, immediate.),
    .valid_FE_rowRanges(x@assays, x@rowRanges, immediate.)
  ))
}

setValidity2("FlexExperiment", function(x) validFE(x, immediate. = FALSE))

## updateObject ################################################################

#' @export
setMethod(
  f = "updateObject",
  signature = "FlexExperiment",
  definition = function(object, ..., verbose = FALSE) {
    object@assays <- updateObject(object@assays, ..., verbose = verbose)
    object@rowRanges <- updateObject(object@rowRanges, ..., verbose = verbose)
    object@colData <- updateObject(object@colData, ..., verbose = verbose)
    object
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FlexExperiment-methods #######################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Methods for FlexExperiment
#'
#' A set of pre-defined accessor methods to extract information from a
#' \code{\link{FlexExperiment}}.
#'
#' @param x,object A \code{\link{FlexExperiment}} object.
#' @param ... `r .dot_param`
#' @param value `r .val_param`
#'
#' @name FlexExperiment-methods
NULL

## show ########################################################################

showFE <- function(object) {
  showS4Title(object)

  cat("---------- assays ----------\n")
  show(object@assays)
  cat("----------------------------\n")

  coolcat("metadata(%d): %s\n", names(metadata(object)))
  coolcat("rowData names(%d): %s\n", names(rowData(object, use.names = FALSE)))
  coolcat("colData names(%d): %s\n", names(colData(object)))
  invisible(NULL)
}

#' @export
#' @rdname FlexExperiment-methods
setMethod("show", "FlexExperiment", showFE)

## length ######################################################################

#' @returns
#' \itemize{
#' \item `length`: A single integer specifying the number of assays.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("length", "FlexExperiment", function(x) length(x@assays))

## names #######################################################################

#' @returns
#' \itemize{
#' \item `names`: A character vector specifying the assay names.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("names", "FlexExperiment", function(x) names(x@assays))

#' @returns
#' \itemize{
#' \item `names<-`: An updated `x` with new assay names.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("names<-", "FlexExperiment", function(x, value) {
  names(x@assays) <- value
  x
})

## dim #########################################################################

#' @returns
#' \itemize{
#' \item `dim`: A two-element vector specifying the row and column numbers of
#' `x`.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("dim", "FlexExperiment", function(x) dim(x@assays))

## dimnames ####################################################################

#' @returns
#' \itemize{
#' \item `dimnames`: A list containing the row and column names of `x`.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("dimnames", "FlexExperiment", function(x) dimnames(x@assays))

#' @returns
#' \itemize{
#' \item `rownames<-` and `colnames<-`: An updated `x` with new dimension names.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("rownames<-", "FlexExperiment", function(x, value) {
  rownames(x@assays) <- value
  x
})

#' @export
#' @rdname FlexExperiment-methods
setMethod("colnames<-", "FlexExperiment", function(x, value) {
  colnames(x@assays) <- value
  x
})

#' @export
#' @rdname FlexExperiment-methods
setMethod("dimnames<-", "FlexExperiment", function(x, value) {
  dimnames(x@assays) <- value
  x
})

## flexible ####################################################################

setMethod("rowFlex", "FlexExperiment", function(x) rowFlex(x@assays))

setMethod("colFlex", "FlexExperiment", function(x) colFlex(x@assays))

setMethod("assayClasses", "FlexExperiment", function(x) assayClasses(x@assays))

## assays ######################################################################

#' @export
#' @rdname FlexExperiment-methods
setMethod("assays", "FlexExperiment", function(x, withDimnames = TRUE, ...) {
  getAssays(x@assays, withDimnames = withDimnames, ...)
})

## assays<- ####################################################################

.always_withDimnames <- function(withDimnames = TRUE) {
  if (withDimnames) {
    return(invisible(NULL))
  }
  msg <- "`withDimnames = FALSE` is always ignored."
  warning(msg, immediate. = TRUE, call. = FALSE)
  invisible(NULL)
}

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assays<-",
  signature = c("FlexExperiment", "FlexAssays"),
  definition = function(x, withDimnames = TRUE, ..., value) {
    .always_withDimnames(withDimnames = withDimnames)
    .set_assays_internal(x, value, ...)
  }
)

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assays<-",
  signature = c("FlexExperiment", "SimpleList"),
  definition = function(x, withDimnames = TRUE, ..., value) {
    .always_withDimnames(withDimnames = withDimnames)
    .set_assays_internal(x, value, ...)
  }
)

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assays<-",
  signature = c("FlexExperiment", "list"),
  definition = function(x, withDimnames = TRUE, ..., value) {
    .always_withDimnames(withDimnames = withDimnames)
    .set_assays_internal(x, value, ...)
  }
)

.set_assays_internal <- function(x, value, ...) {
  value <- FlexAssays(
    value,
    rowFlex = rowFlex(x),
    colFlex = colFlex(x),
    assayClasses = assayClasses(x)
  )
  x@assays <- value
  validFE(x)
  x
}

## assay #######################################################################

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assay",
  signature = c("FlexExperiment", "missing"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    if (length(x) == 0) {
      fmt <- "'assay(<%s>, i=\"missing\", ...)' failed: no avaliable assay"
      stop(sprintf(fmt, class(x)[1]))
    }
    getOneAssay(x@assays, 1L, withDimnames = withDimnames)
  }
)

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assay",
  signature = c("FlexExperiment", "integer"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    .get_assay_internal(x, i, withDimnames = withDimnames, ...)
  }
)

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assay",
  signature = c("FlexExperiment", "numeric"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    .get_assay_internal(x, i, withDimnames = withDimnames, ...)
  }
)

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assay",
  signature = c("FlexExperiment", "character"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    .get_assay_internal(x, i, withDimnames = withDimnames, ...)
  }
)

.get_assay_internal <- function(x, i, withDimnames = TRUE, ...) {
  tryCatch(
    expr = getOneAssay(x@assays, i, withDimnames = withDimnames),
    error = function(err) {
      e <- "'assay(<%s>, i = <%s>, ...)' failed. Invalid subscript 'i': %s\n"
      stop(sprintf(e, class(x)[1], mode(i)[1], i), conditionMessage(err))
    }
  )
}

## assay<- #####################################################################

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assay<-",
  signature = c("FlexExperiment", "missing"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    i <- 1L
    .always_withDimnames(withDimnames = withDimnames)
    .set_assay_internal(x = x, i = i, value = value, ...)
  }
)

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assay<-",
  signature = c("FlexExperiment", "integer"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    .always_withDimnames(withDimnames = withDimnames)
    .set_assay_internal(x = x, i = i, value = value, ...)
  }
)

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assay<-",
  signature = c("FlexExperiment", "numeric"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    .always_withDimnames(withDimnames = withDimnames)
    .set_assay_internal(x = x, i = i, value = value, ...)
  }
)

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "assay<-",
  signature = c("FlexExperiment", "character"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    .always_withDimnames(withDimnames = withDimnames)
    .set_assay_internal(x = x, i = i, value = value, ...)
  }
)

.set_assay_internal <- function(x, i, value, ...) {
  x@assays <- setOneAssay(x@assays, i, value)
  return(x)
}

## assayNames ##################################################################

#' @export
#' @rdname FlexExperiment-methods
setMethod("assayNames", "FlexExperiment", function(x) names(x@assays))

## assayNames<- ################################################################

#' @export
#' @rdname FlexExperiment-methods
setMethod("assayNames<-", "FlexExperiment", function(x, value) {
  names(x@assays) <- value
  x
})

## rowRanges ###################################################################

#' @export
#' @rdname FlexExperiment-methods
setMethod("rowRanges", "FlexExperiment", function(x, ...) {
  rr <- updateObject(x@rowRanges, check = FALSE)
  names(rr) <- rownames(x)
  rr
})

## rowRanges<- #################################################################

#' @export
#' @rdname FlexExperiment-methods
setMethod("rowRanges<-", c("FlexExperiment", "NULL"), function(x, ..., value) {
  fmt <- "Cannot set rowRanges(<%s>) to NULL."
  warning(sprintf(fmt, class(x)[1]), immediate. = TRUE, call. = FALSE)
  x
})

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "rowRanges<-",
  signature = c("FlexExperiment", "GenomicRanges"),
  definition = function(x, ..., value) .set_rowRanges_FE(x, value)
)

#' @export
#' @rdname FlexExperiment-methods
setMethod(
  f = "rowRanges<-",
  signature = c("FlexExperiment", "GRangesList"),
  definition = function(x, ..., value) .set_rowRanges_FE(x, value)
)

.set_rowRanges_FE <- function(x, value) {
  value <- value[rownames(x)]
  .valid_FE_rowRanges(x@assays, value)
  names(value) <- NULL
  slot(x, "rowRanges", check = FALSE) <- value
  x
}

## rowData #####################################################################

#' @param use.names Whether or not to set the dimension names of `x` as the row
#' names of output `r .doc_links("DataFrame")`.
#'
#' @returns
#' \itemize{
#' \item `rowData`: A `r .doc_links("DataFrame")` specifying the meta data for
#' rows of `x`.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("rowData", "FlexExperiment", function(x, use.names = TRUE, ...) {
  df <- mcols(x@rowRanges, use.names = FALSE, ...)
  if (use.names) {
    rownames(df) <- rownames(x)
  }
  df
})

## rowData<- ###################################################################

#' @returns
#' \itemize{
#' \item `rowData<-`: An updated `x` with new metadata of rows.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("rowData<-", "FlexExperiment", function(x, ..., value) {
  value <- .prep_dimData(x@assays, value, nrow, rownames)
  r <- rowRanges(x)
  mcols(r) <- value
  slot(x, "rowRanges", check = FALSE) <- r
  x
})

.prep_dimData <- function(assays, dimData, nfun, namefun) {
  if (is.null(dimData)) {
    return(new2("DFrame", nrows = nfun(assays), check = FALSE))
  }
  if (!inherits(dimData, "DataFrame")) {
    dimData <- as(dimData, "DataFrame")
  }
  err <- sprintf("<%s> row names out of bound: ", class(dimData)[1])
  idx <- getValidCharBound(namefun(assays), rownames(dimData), err = err)
  if (length(idx) != nfun(assays)) {
    funstr <- as.character(substitute(namefun))
    stop(sprintf("rownames(value) does not match %s(x).", funstr))
  }
  dimData <- updateObject(dimData, check = FALSE)[idx, , drop = FALSE]
  rownames(dimData) <- NULL
  dimData
}

## colData #####################################################################

#' @returns
#' \itemize{
#' \item `colData`: A `r .doc_links("DataFrame")` specifying the meta data for
#' columns of `x`.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("colData", "FlexExperiment", function(x, ...) {
  df <- updateObject(x@colData, check = FALSE)
  rownames(df) <- colnames(x)
  df
})

## colData<- ###################################################################

#' @returns
#' \itemize{
#' \item `colData<-`: An updated `x` with new metadata of columns.
#' }
#'
#' @export
#' @rdname FlexExperiment-methods
setMethod("colData<-", "FlexExperiment", function(x, ..., value) {
  value <- .prep_dimData(x@assays, value, ncol, colnames)
  slot(x, "colData", check = FALSE) <- value
  x
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FlexExperiment-subset ########################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Extraction methods for FlexExperiment
#'
#' Methods for extraction or subsetting \code{\link{FlexExperiment}} with
#' `[`, `[[` or `$`.
#'
#' @param x A \code{\link{FlexExperiment}} object.
#' @param ... `r .dot_param`
#' @param value `r .val_param`
#'
#' @name FlexExperiment-subset
NULL

## `[` #########################################################################

#' @param i,j Indices specifying elements to extract or replace.
#' @param drop Whether or not to drop the empty assays when subset `x`.
#'
#' @returns
#' \itemize{
#' \item `[`: A subset `x`.
#' }
#'
#' @export
#' @rdname FlexExperiment-subset
setMethod("[", "FlexExperiment", function(x, i, j, ..., drop = FALSE) {
  miss.i <- missing(i)
  miss.j <- missing(j)
  if (miss.i & miss.j) {
    return(x)
  }
  if (!miss.i) {
    ii <- .subset_names2idx(i, rownames(x))
    x@rowRanges <- x@rowRanges[ii]
    if (miss.j) {
      x@assays <- x@assays[i, , drop = drop]
      return(x)
    }
  }
  if (!miss.j) {
    jj <- .subset_names2idx(j, colnames(x))
    x@colData <- x@colData[jj, , drop = FALSE]
    if (miss.i) {
      x@assays <- x@assays[, j, drop = drop]
      return(x)
    }
  }
  x@assays <- x@assays[i, j, drop = drop]
  return(x)
})

.subset_names2idx <- function(subset, names) {
  if (!is.character(subset)) {
    return(as.vector(subset))
  }
  getValidCharBound(subset, names, err = "out of bound: ")
}

## `[[` ########################################################################

#' @returns
#' \itemize{
#' \item `[[` or `$`: One column of `colData(x)`.
#' }
#'
#' @export
#' @rdname FlexExperiment-subset
setMethod("[[", "FlexExperiment", function(x, i, j, ...) colData(x)[[i, ...]])

## `[[<-` ######################################################################

#' @returns
#' \itemize{
#' \item `[[<-` or `$<-`: An updated `x` with `value` added to one specific
#' column of `colData(x)`.
#' }
#'
#' @export
#' @rdname FlexExperiment-subset
setMethod("[[<-", "FlexExperiment", function(x, i, j, ..., value) {
  colData(x)[[i, ...]] <- value
  x
})

## .DollarNames ################################################################

#' @param patten A regular expression. Only matching names are returned.
#'
#' @method .DollarNames FlexExperiment
#' @export
#' @rdname FlexExperiment-subset
.DollarNames.FlexExperiment <- function(x, pattern = "") {
  grep(pattern, names(colData(x)), value = TRUE)
}

#' @method $ FlexExperiment
#' @export
#' @rdname FlexExperiment-subset
"$.FlexExperiment" <- function(x, name) colData(x)[[name]]

#' @method $<- FlexExperiment
#' @export
#' @rdname FlexExperiment-subset
"$<-.FlexExperiment" <- function(x, name, value) {
  colData(x)[[name]] <- value
  x
}

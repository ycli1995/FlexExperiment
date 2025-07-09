#' @include FlexExperiment-class.R
#' @include SingleCellFlexExperiment-class.R
#'
#' @importFrom SingleCellExperiment altExp altExp<- altExpNames altExpNames<-
#' altExps<-
#' @importFrom rlang is_bare_list
#' @importFrom methods as
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom Matrix dgCMatrix
NULL

validateConcatData <- function(obj.list, namefun) {
  for (i in seq_along(obj.list)) {
    if (anyDuplicated(namefun(obj.list[[i]]))) {
      funstr <- as.character(substitute(namefun))
      stop("The input object list [[", i, "]] contains duplicated ", funstr)
    }
  }
  invisible(NULL)
}

#' Concatenate matrices
#'
#' Methods to concatenate matrix-like objects by rows or columns.
#'
#' @param x A matrix-like object
#' @param y A single or a list of matrix-like object(s) to be merged.
#' @param ... `r .dot_param`
#'
#' @details
#' When concatenate matrices by rows, the union of original column names will be
#' used. Any duplicated row names among matrices will be enforced unique by
#' adding suffixes with \code{\link{make.unique}}.
#'
#' When concatenate matrices by columns, use the opposite way.
#'
#' @returns
#' A concatenated matrix-like object.
#'
#' @name concat-matrices
NULL

#' @rdname concat-matrices
#' @export concatByRows
concatByRows <- function(x, y, ...) UseMethod("concatByRows", x)

#' @rdname concat-matrices
#' @export concatByCols
concatByCols <- function(x, y, ...) UseMethod("concatByCols", x)

#' @rdname concat-matrices
#' @export
#' @method concatByRows matrix
concatByRows.matrix <- function(x, y, ...) {
  if (!is_bare_list(y)) {
    y <- list(y)
  }
  mat.list <- c(list(x), y)
  validateConcatData(mat.list, colnames)
  .concat_matrix(mat.list = mat.list, by = "row")
}

#' @rdname concat-matrices
#' @export
#' @method concatByCols matrix
concatByCols.matrix <- function(x, y, ...) {
  if (!is_bare_list(y)) {
    y <- list(y)
  }
  mat.list <- c(list(x), y)
  validateConcatData(mat.list, rownames)
  .concat_matrix(mat.list = mat.list, by = "column")
}

#' @rdname concat-matrices
#' @export
#' @method concatByRows dgCMatrix
concatByRows.dgCMatrix <- function(x, y, ...) {
  if (!is_bare_list(y)) {
    y <- list(y)
  }
  mat.list <- c(list(x), y)
  validateConcatData(mat.list, colnames)
  .concat_dgcmatrix(mat.list = mat.list, by = "row")
}

#' @rdname concat-matrices
#' @export
#' @method concatByCols dgCMatrix
concatByCols.dgCMatrix <- function(x, y, ...) {
  if (!is_bare_list(y)) {
    y <- list(y)
  }
  mat.list <- c(list(x), y)
  validateConcatData(mat.list, rownames)
  .concat_dgcmatrix(mat.list = mat.list, by = "column")
}

#' @export
#' @method concatByRows data.frame
concatByRows.data.frame <- function(x, y, join = c("inner", "outer"), ...) {
  join <- match.arg(join)
  if (!is_bare_list(y)) {
    y <- list(y)
  }
  df.list <- c(list(x), y)
  .concat_DFs_by_rows(df.list, join = join, ...)
}

#' @method concatByRows DataFrame
concatByRows.DataFrame <- function(x, y, join = c("inner", "outer"), ...) {
  join <- match.arg(join)
  concatByRows.data.frame(x = x, y = y, join = join, ...)
}

#' Concatenate complex matrix-like objects
#'
#' Methods to concatenate matrix-like objects to merge omics data.
#'
#' @inheritParams concat-matrices
#'
#' @returns
#' A concatenated matrix-like object.
#'
#' @name concatByCols-objects
NULL

#' @param assays Names of `r .doc_links("assays")` to be merged. If set to
#' `NULL`, common assays will be automatically searched from input objects. If
#' set to `FALSE`, no assay will be merged.
#' @param label Name of a new column in `colData()` to hold the batch
#' information. Default is `NULL` (No batch information added).
#' @param add.idx A character vector with `length(c(x, y))` elements; appends
#' the corresponding values to the start of each objects' cell names.
#' @param collapse String to collapse `add.idx` and cell names of each input
#' object. Default is `"_"`.
#' @param verbose `r .vb_param`
#'
#' @rdname concatByCols-objects
#' @export
#' @method concatByCols SummarizedExperiment
concatByCols.SummarizedExperiment <- function(
    x, y,
    assays = NULL,
    label = NULL,
    add.idx = NULL,
    collapse = "_",
    verbose = TRUE,
    ...
) {
  SCEs <- .prep_merge_xy(x = x, y = y, class = "SummarizedExperiment")
  if (length(SCEs) == 1) {
    return(SCEs[[1]])
  }
  SCEs <- .prep_merge_SCEs(
    SCEs,
    label = label,
    add.idx = add.idx,
    collapse = collapse,
    verbose = verbose
  )
  # colData and rowData
  verboseMsg("-----> Merging rowData...")
  rdata <- .concat_DFs(lapply(SCEs, rowData), "outer", merge.by = "first")
  verboseMsg("-----> Merging colData...")
  cdata <- .concat_DFs_by_rows(lapply(SCEs, colData), join = "outer")

  # merge assays
  new.assays <- .merge_SCEs_assays(
    SCEs = SCEs,
    assays = assays,
    verbose = verbose
  )
  SummarizedExperiment(
    assays = new.assays,
    rowData = rdata,
    colData = cdata
  )
}

#' @param altExps Names of `r .doc_links("altExps")` to be merged. If set to
#' `NULL`, common alternative experiments will be automatically searched from
#' input objects. If set to `FALSE`, no `altExps` will be merged.
#' @param altExps.assays Names of `assays` to be merged in `altExps`.
#'
#' @rdname concatByCols-objects
#' @export
#' @method concatByCols SingleCellExperiment
concatByCols.SingleCellExperiment <- function(
    x, y,
    assays = NULL,
    reducedDims = NULL,
    altExps = NULL,
    altExps.assays = NULL,
    label = NULL,
    add.idx = NULL,
    collapse = "_",
    verbose = TRUE,
    ...
) {
  SCEs <- .prep_merge_xy(x = x, y = y, class = "SingleCellExperiment")
  if (length(SCEs) == 1) {
    return(SCEs[[1]])
  }
  SCEs <- .prep_merge_SCEs(
    SCEs,
    label = label,
    add.idx = add.idx,
    collapse = collapse,
    verbose = verbose
  )
  # colData and rowData
  verboseMsg("-----> Merging rowData...")
  rdata <- .concat_DFs(lapply(SCEs, rowData), "outer", merge.by = "first")
  verboseMsg("-----> Merging colData...")
  cdata <- .concat_DFs_by_rows(lapply(SCEs, colData), join = "outer")

  # merge altExps
  new.altExps <- .merge_SCEs_altExps(
    SCEs = SCEs,
    altExps = altExps,
    altExps.assays = altExps.assays,
    verbose = verbose
  )
  # merge assays
  new.assays <- .merge_SCEs_assays(
    SCEs = SCEs,
    assays = assays,
    verbose = verbose
  )
  # merge reductions
  new.reducs <- .merge_SCEs_reducedDims(
    SCEs = SCEs,
    reducedDims = reducedDims,
    verbose = verbose
  )
  new.SCE <- SingleCellExperiment(
    assays = new.assays,
    rowData = rdata,
    colData = cdata,
    reducedDims = new.reducs,
    altExps = new.altExps
  )
  new.SCE
}

#' @param assayClasses,rowFlex,colFlex Parameters passed to the new merged
#' \code{\link{FlexExperiment}}.
#'
#' @rdname concatByCols-objects
#' @export
#' @method concatByCols FlexExperiment
concatByCols.FlexExperiment <- function(
    x, y,
    assays = NULL,
    assayClasses = NULL,
    rowFlex = "bounded",
    colFlex = "bounded",
    label = NULL,
    add.idx = NULL,
    collapse = "_",
    verbose = TRUE,
    ...
) {
  SCEs <- .prep_merge_xy(x = x, y = y, class = "FlexExperiment")
  if (length(SCEs) == 1) {
    return(SCEs[[1]])
  }
  SCEs <- .prep_merge_SCEs(
    SCEs,
    label = label,
    add.idx = add.idx,
    collapse = collapse,
    verbose = verbose
  )

  # colData and rowData
  verboseMsg("-----> Merging rowData...")
  rdata <- .concat_DFs(lapply(SCEs, rowData), "outer", merge.by = "first")
  verboseMsg("-----> Merging colData...")
  cdata <- .concat_DFs_by_rows(lapply(SCEs, colData), join = "outer")
  new.assays <- .merge_SCEs_assays(
    SCEs = SCEs,
    assays = assays,
    verbose = verbose
  )
  FlexExperiment(
    assays = new.assays,
    rowData = rdata,
    colData = cdata,
    assayClasses = assayClasses,
    rowFlex = rowFlex,
    colFlex = colFlex
  )
}

#' @rdname concatByCols-objects
#' @export
#' @method concatByCols SingleCellFlexExperiment
concatByCols.SingleCellFlexExperiment <- function(
    x, y,
    assays = NULL,
    reducedDims = NULL,
    assayClasses = NULL,
    rowFlex = "bounded",
    colFlex = "bounded",
    label = NULL,
    add.idx = NULL,
    collapse = "_",
    verbose = TRUE,
    ...
) {
  SCEs <- .prep_merge_xy(x = x, y = y, class = "FlexExperiment")
  if (length(SCEs) == 1) {
    return(SCEs[[1]])
  }
  SCEs <- .prep_merge_SCEs(
    SCEs,
    label = label,
    add.idx = add.idx,
    collapse = collapse,
    verbose = verbose
  )

  # colData and rowData
  verboseMsg("-----> Merging rowData...")
  rdata <- .concat_DFs(lapply(SCEs, rowData), "outer", merge.by = "first")
  verboseMsg("-----> Merging colData...")
  cdata <- .concat_DFs_by_rows(lapply(SCEs, colData), join = "outer")
  new.assays <- .merge_SCEs_assays(
    SCEs = SCEs,
    assays = assays,
    verbose = verbose
  )

  # merge reductions
  new.reducs <- .merge_SCEs_reducedDims(
    SCEs = SCEs,
    reducedDims = reducedDims,
    verbose = verbose
  )
  SingleCellFlexExperiment(
    assays = new.assays,
    rowData = rdata,
    colData = cdata,
    reducedDims = new.reducs,
    assayClasses = assayClasses,
    rowFlex = rowFlex,
    colFlex = colFlex
  )
}

.check_same <- function(obj.list) {
  vapply(
    X = obj.list,
    FUN = all.equal,
    FUN.VALUE = logical(1L),
    current = list[[1]],
    check.attributes = FALSE
  )
}

.check_null <- function(list) {
  vapply(list, is.null, FUN.VALUE = logical(1L))
}

.get_valid_merge_by <- function(merge.by = NULL) {
  if (length(merge.by) == 0) {
    return(merge.by)
  }
  if (length(merge.by) > 1) {
    warning(
      "Only take the first 'merge.by': ", merge.by[1],
      immediate. = TRUE, call. = FALSE
    )
    merge.by <- merge.by[1]
  }
  .merge_by <- c("first")
  if (merge.by %in% .merge_by) {
    return(merge.by)
  }
  .merge_by <- paste(.merge_by, collapse = ", ")
  stop("'merge.by' should only be one of the following: \n  ", .merge_by)
}

.concat_DFs_helper <- function(df.list, new.df, merge.by) {
  merge.by <- .get_valid_merge_by(merge.by = merge.by)
  if (is.null(merge.by)) {
    return(new.df)
  }
  use.rows <- rownames(new.df)
  use.cols <- unique(unlist(lapply(df.list, colnames)))
  if (merge.by == "first") {
    # The first column seen at each from each position.
    if (length(use.cols) == 0) {
      return(new.df)
    }
    new.df[, use.cols] <- NA
    for (i in rev(seq_along(df.list))) {
      tmp.rows <- intersect(use.rows, rownames(df.list[[i]]))
      tmp.cols <- intersect(use.cols, colnames(df.list[[i]]))
      new.df[tmp.rows, tmp.cols] <- df.list[[i]][tmp.rows, tmp.cols]
    }
    return(new.df)
  }
  return(new.df)
}

.concat_DFs <- function(
    df.list,
    r.join = c("inner", "outer"),
    merge.by = NULL
) {
  r.join <- match.arg(r.join)
  all.rows <- lapply(df.list, rownames)
  all.cols <- lapply(df.list, colnames)
  if (r.join == "inner") {
    use.rows <- Reduce(intersect, all.rows)
    new.df <- df.list[[1]][use.rows, character(), drop = FALSE]
  } else {
    all.rows <- unlist(all.rows)
    new.df <- Reduce(
      f = rbind,
      x = lapply(df.list, function(df) df[, character(), drop = FALSE])
    )
    new.df <- new.df[!duplicated(all.rows), , drop = FALSE]
    rownames(new.df) <- unique(all.rows)
  }
  .concat_DFs_helper(df.list, new.df, merge.by)
}

.concat_DFs_by_rows <- function(df.list, join = c("inner", "outer"), ...) {
  join <- match.arg(join)
  all.cols <- lapply(df.list, FUN = colnames)
  if (join == "inner") {
    use.cols <- Reduce(intersect, all.cols)
    for (i in seq_along(df.list)) {
      df.list[[i]] <- df.list[[i]][, use.cols, drop = FALSE]
    }
    new.df <- Reduce(rbind, df.list)
    return(new.df)
  }
  use.cols <- unique(unlist(all.cols))
  for (i in seq_along(df.list)) {
    na.cols <- setdiff(use.cols, colnames(df.list[[i]]))
    if (length(na.cols) > 0) {
      df.list[[i]][, na.cols] <- NA
    }
  }
  new.df <- Reduce(rbind, df.list)
  new.df
}


.concat_matrix <- function(mat.list, by = c("row", "column")) {
  by <- match.arg(by)
  all.colnames <- lapply(mat.list, colnames)
  all.rownames <- lapply(mat.list, rownames)
  if (by == "row") {
    uniq.colnames <- unique(unlist(all.colnames))
    uniq.rownames <- make.unique(unlist(all.rownames))
  } else {
    uniq.colnames <- make.unique(unlist(all.colnames))
    uniq.rownames <- unique(unlist(all.rownames))
  }
  m <- matrix(0, nrow = length(uniq.rownames), ncol = length(uniq.colnames))
  rownames(m) <- uniq.rownames
  colnames(m) <- uniq.colnames
  i1 <- 0
  if (by == "row") {
    for (i in seq_along(mat.list)) {
      i1 <- i1 + 1
      i2 <- i1 + nrow(mat.list[[i]]) - 1
      m[i1:i2, all.colnames[[i]]] <- as.matrix(mat.list[[i]])
      i1 <- i2
    }
  } else {
    for (i in seq_along(mat.list)) {
      i1 <- i1 + 1
      i2 <- i1 + ncol(mat.list[[i]]) - 1
      m[all.rownames[[i]], i1:i2] <- as.matrix(mat.list[[i]])
      i1 <- i2
    }
  }
  return(m)
}

.concat_dgcmatrix <- function(mat.list, by = c("row", "column")) {
  by <- match.arg(by)
  for (i in seq_along(mat.list)) {
    if (!inherits(mat.list[[i]], "dgCMatrix")) {
      mat.list[[i]] <- as(mat.list[[i]], "dgCMatrix")
    }
  }
  all.colnames <- lapply(mat.list, colnames)
  all.rownames <- lapply(mat.list, rownames)
  if (by == "row") {
    uniq.colnames <- unique(unlist(all.colnames))
    uniq.rownames <- make.unique(unlist(all.rownames))
    use.rbind <- all(duplicated(all.colnames)[-1])
    if (use.rbind) {
      new.mat <- Reduce(rbind, mat.list)
    } else {
      new.mat <- row_merge_dgcmatrix_cpp(mat.list, all.colnames, uniq.colnames)
      colnames(new.mat) <- uniq.colnames
    }
    rownames(new.mat) <- uniq.rownames
  } else {
    uniq.colnames <- make.unique(unlist(all.colnames))
    uniq.rownames <- unique(unlist(all.rownames))
    use.cbind <- all(duplicated(all.rownames)[-1])
    if (use.cbind) {
      new.mat <- Reduce(cbind, mat.list)
    } else {
      new.mat <- col_merge_dgcmatrix_cpp(mat.list, all.rownames, uniq.rownames)
      rownames(new.mat) <- uniq.rownames
    }
    colnames(new.mat) <- uniq.colnames
  }
  new.mat
}

.prep_merge_xy <- function(x, y, class) {
  if (!is_bare_list(y)) {
    y <- list(y)
  }
  SCEs <- c(list(x), y)
  checkInherits(SCEs, class = class)
  SCEs
}

.prep_merge_SCEs <- function(
    SCEs,
    label = NULL,
    add.idx = NULL,
    collapse = "_",
    verbose = TRUE
) {
  verboseMsg("-----> Preparing ", class(SCEs[[1]])[1], "(s) to be merged...")
  if (!is.null(add.idx)) {
    stopifnot(length(add.idx) == length(SCEs))
    verboseMsg("Adding idx prefixes: ", paste(add.idx, collapse = ", "))
    for (i in seq_along(SCEs)) {
      colnames(SCEs[[i]]) <- paste0(add.idx[i], collapse, colnames(SCEs[[i]]))
    }
  }
  if (!is.null(label)) {
    label <- label[1]
    if (length(names(SCEs)) == 0) {
      names(SCEs) <- as.character(seq_along(SCEs))
    }
    verboseMsg(
      "Adding column '", label, "' to specify batch info: ",
      paste(names(SCEs), collapse = ", ")
    )
    for (i in names(SCEs)) {
      SCEs[[i]][[label]] <- i
    }
  }
  SCEs <- dedupNames(SCEs, getfunc = colnames, setfunc = `colnames<-`)
  return(SCEs)
}

.merge_SCEs_reducedDims <- function(
    SCEs,
    reducedDims = NULL,
    verbose = TRUE
) {
  verboseMsg("-----> Merging reducedDims...")
  reducs <- searchCommonNames(
    obj.list = SCEs,
    names.func = reducedDimNames,
    names = reducedDims,
    verbose = verbose
  )
  new.reds <- list()
  for (i in reducs) {
    verboseMsg("Merging '", i, "'")
    old.reds <- list()
    for (j in seq_along(SCEs)) {
      old.reds[[j]] <- as.matrix(reducedDim(SCEs[[j]], i, withDimnames = TRUE))
    }
    min.ncol <- vapply(
      X = old.reds,
      FUN = ncol,
      FUN.VALUE = integer(1L),
      USE.NAMES = FALSE
    )
    if (length(unique(min.ncol)) > 1) {
      min.ncol <- min(min.ncol)
      warning(
        "reducedDim '", i, "' contain differing numbers of dimensions, ",
        "merging the first ", min.ncol,
        call. = FALSE, immediate. = TRUE
      )
      for (j in seq_along(old.reds)) {
        old.reds[[j]] <- old.reds[[j]][, seq_len(min.ncol)]
      }
    }
    new.reds[[i]] <- as(Reduce(rbind, old.reds), "LinearEmbeddingMatrix")
  }
  new.reds
}

.merge_SCEs_assays <- function(SCEs, assays = NULL, verbose = TRUE) {
  verboseMsg("-----> Merging assays...")
  assays <- searchCommonNames(
    obj.list = SCEs,
    names.func = assayNames,
    names = assays,
    verbose = verbose
  )
  new.assays <- list()
  for (i in assays) {
    verboseMsg("Merging '", i, "'")
    old.assays <- list()
    for (j in seq_along(SCEs)) {
      old.assays[[j]] <- assay(SCEs[[j]], i = i, withDimnames = TRUE)
    }
    new.assays[[i]] <- concatByCols(old.assays[[1]], old.assays[-1])
  }
  new.assays
}

.merge_SCEs_altExps <- function(
    SCEs,
    altExps = NULL,
    altExps.assays = NULL,
    verbose = TRUE
) {
  verboseMsg("-----> Merging altExp...")
  altExps <- searchCommonNames(
    obj.list = SCEs,
    names.func = altExpNames,
    names = altExps,
    verbose = verbose
  )
  new.altExps <- list()
  for (i in altExps) {
    verboseMsg("Merging '", i, "'")
    old.altExps <- list()
    for (j in seq_along(SCEs)) {
      old.altExps[[j]] <- altExp(SCEs[[j]], e = i, withDimnames = TRUE)
      colData(old.altExps[[j]]) <- NULL
      altExps(old.altExps[[j]]) <- NULL
      reducedDims(old.altExps[[j]]) <- NULL
    }
    new.altExps[[i]] <- concatByCols(
      x = old.altExps[[1]],
      y = old.altExps[2:length(old.altExps)],
      assays = altExps.assays,
      reducedDims = FALSE,
      altExps = FALSE,
      label = NULL,
      add.idx = NULL,
      verbose = verbose
    )
  }
  new.altExps
}





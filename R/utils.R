#' @importFrom FlexAssays getErrors identicalNoAttrs
NULL

.identical_fmatch <- function(x, y, ...) {
  identicalNoAttrs(x, y, ignore.attrs = ".match.hash", ...)
}

#' Validate Character Vector Membership
#'
#' Checks whether all elements of a character vector are present in a reference
#' vector, and optionally returns their index positions. Throws an error if any
#' element is not found.
#'
#' @param query A character vector. The values to be matched against `ref`.
#' @param ref A character vector. The reference set of allowed values.
#' @param return Logical (default `TRUE`). If `TRUE`, return the index positions
#' of `query` in `ref`; if `FALSE`, return `NULL` invisibly.
#' @param err A character string. An optional prefix message to include in the
#' error if `query` has invalid elements.
#'
#' @return If `return = TRUE`, an integer vector of indices indicating the
#' position of each element of `query` in `ref`. If all values are matched but
#' `return = FALSE`, returns `invisible(NULL)`.
#'
#' @examples
#' getValidCharBound(c("a", "b"), letters)      # returns c(1, 2)
#' getValidCharBound(c("x", "z"), letters)      # returns c(24, 26)
#' # getValidCharBound(c("a", "A"), letters)    # throws error
#'
#' @export
getValidCharBound <- function(query, ref, return = TRUE, err = "") {
  idx <- match(query, ref)
  if (any(bad <- is.na(idx))) {
    stop(err, "\n ", paste(query[head(which(bad))], collapse = ", "), "...")
  }
  if (return) {
    return(idx)
  }
  return(invisible(NULL))
}

#' Check Identical Objects with Custom Error Handling
#'
#' Compares two R objects for identity using a custom comparison function. If
#' the objects are not identical, an error is raised using a user-specified
#' message.
#'
#' @param x,y Objects to compare.
#' @param err A character string. Error message to show if `x` and `y` are not
#' identical.
#' @param immediate. Logical (default `TRUE`). Whether to immediately signal the
#' error (passed to `getErrors()`).
#'
#' @return Returns `invisible(NULL)` if `x` and `y` are identical. Otherwise, an
#' error is thrown.
#'
#' @details This function is a wrapper for equality checking that allows for
#' custom error messages and deferred error signaling. It uses `getErrors()` to
#' delegate error generation.
#'
#' @examples
#' # checkIdentical(1:5, 1:5)  # passes silently
#' # checkIdentical(1:5, 1:4)  # throws error
#'
#' @export
checkIdentical <- function(x, y, err = "", immediate. = TRUE) {
  if (.identical_fmatch(x, y)) {
    return(invisible(NULL))
  }
  getErrors(err, immediate. = immediate.)
}

#' Check Equality of Scalar Values with Custom Error
#'
#' Compares two scalar values for equality. If they are not equal, throws a
#' custom error.
#'
#' @param x,y Scalar values to compare.
#' @param err A character string. Custom error message to show if `x` and `y`
#' are not equal.
#' @param immediate. Logical (default `TRUE`). Whether to immediately signal
#' the error (passed to `getErrors()`).
#'
#' @return Returns `invisible(NULL)` if `x == y`. Otherwise, throws an error
#' with the given message.
#'
#' @details This function is intended for use in internal validation checks.
#'
#' @examples
#' # checkEqual(1, 1)     # passes silently
#' # checkEqual(1, 2)     # throws error
#'
#' @keywords internal
checkEqual <- function(x, y, err = "", immediate. = TRUE) {
  if (x == y) {
    return(invisible(NULL))
  }
  getErrors(err, immediate. = immediate.)
}

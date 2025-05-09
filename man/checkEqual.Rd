% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{checkEqual}
\alias{checkEqual}
\title{Check Equality of Scalar Values with Custom Error}
\usage{
checkEqual(x, y, err = "", immediate. = TRUE)
}
\arguments{
\item{x, y}{Scalar values to compare.}

\item{err}{A character string. Custom error message to show if \code{x} and \code{y}
are not equal.}

\item{immediate.}{Logical (default \code{TRUE}). Whether to immediately signal
the error (passed to \code{getErrors()}).}
}
\value{
Returns \code{invisible(NULL)} if \code{x == y}. Otherwise, throws an error
with the given message.
}
\description{
Compares two scalar values for equality. If they are not equal, throws a
custom error.
}
\details{
This function is intended for use in internal validation checks.
}
\examples{
# checkEqual(1, 1)     # passes silently
# checkEqual(1, 2)     # throws error

}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FlexExperiment-class.R
\name{FlexExperiment-methods}
\alias{FlexExperiment-methods}
\alias{show,FlexExperiment-method}
\alias{length,FlexExperiment-method}
\alias{names,FlexExperiment-method}
\alias{names<-,FlexExperiment,ANY-method}
\alias{dim,FlexExperiment-method}
\alias{dimnames,FlexExperiment-method}
\alias{rownames<-,FlexExperiment-method}
\alias{colnames<-,FlexExperiment-method}
\alias{dimnames<-,FlexExperiment,ANY-method}
\alias{assays,FlexExperiment-method}
\alias{assays<-,FlexExperiment,FlexAssays-method}
\alias{assays<-,FlexExperiment,SimpleList-method}
\alias{assays<-,FlexExperiment,list-method}
\alias{assay,FlexExperiment,missing-method}
\alias{assay,FlexExperiment,integer-method}
\alias{assay,FlexExperiment,numeric-method}
\alias{assay,FlexExperiment,character-method}
\alias{assay<-,FlexExperiment,missing-method}
\alias{assay<-,FlexExperiment,integer-method}
\alias{assay<-,FlexExperiment,numeric-method}
\alias{assay<-,FlexExperiment,character-method}
\alias{assayNames,FlexExperiment-method}
\alias{assayNames<-,FlexExperiment,ANY-method}
\alias{rowRanges,FlexExperiment-method}
\alias{rowRanges<-,FlexExperiment,NULL-method}
\alias{rowRanges<-,FlexExperiment,GenomicRanges-method}
\alias{rowRanges<-,FlexExperiment,GRangesList-method}
\alias{rowData,FlexExperiment-method}
\alias{rowData<-,FlexExperiment-method}
\alias{colData,FlexExperiment-method}
\alias{colData<-,FlexExperiment,ANY-method}
\title{Methods for FlexExperiment}
\usage{
\S4method{show}{FlexExperiment}(object)

\S4method{length}{FlexExperiment}(x)

\S4method{names}{FlexExperiment}(x)

\S4method{names}{FlexExperiment,ANY}(x) <- value

\S4method{dim}{FlexExperiment}(x)

\S4method{dimnames}{FlexExperiment}(x)

\S4method{rownames}{FlexExperiment}(x) <- value

\S4method{colnames}{FlexExperiment}(x) <- value

\S4method{dimnames}{FlexExperiment,ANY}(x) <- value

\S4method{assays}{FlexExperiment}(x, withDimnames = TRUE, ...)

\S4method{assays}{FlexExperiment,FlexAssays}(x, withDimnames = TRUE, ...) <- value

\S4method{assays}{FlexExperiment,SimpleList}(x, withDimnames = TRUE, ...) <- value

\S4method{assays}{FlexExperiment,list}(x, withDimnames = TRUE, ...) <- value

\S4method{assay}{FlexExperiment,missing}(x, i, withDimnames = TRUE, ...)

\S4method{assay}{FlexExperiment,integer}(x, i, withDimnames = TRUE, ...)

\S4method{assay}{FlexExperiment,numeric}(x, i, withDimnames = TRUE, ...)

\S4method{assay}{FlexExperiment,character}(x, i, withDimnames = TRUE, ...)

\S4method{assay}{FlexExperiment,missing}(x, i, withDimnames = TRUE, ...) <- value

\S4method{assay}{FlexExperiment,integer}(x, i, withDimnames = TRUE, ...) <- value

\S4method{assay}{FlexExperiment,numeric}(x, i, withDimnames = TRUE, ...) <- value

\S4method{assay}{FlexExperiment,character}(x, i, withDimnames = TRUE, ...) <- value

\S4method{assayNames}{FlexExperiment}(x)

\S4method{assayNames}{FlexExperiment,ANY}(x) <- value

\S4method{rowRanges}{FlexExperiment}(x, ...)

\S4method{rowRanges}{FlexExperiment,NULL}(x, ...) <- value

\S4method{rowRanges}{FlexExperiment,GenomicRanges}(x, ...) <- value

\S4method{rowRanges}{FlexExperiment,GRangesList}(x, ...) <- value

\S4method{rowData}{FlexExperiment}(x, use.names = TRUE, ...)

\S4method{rowData}{FlexExperiment}(x, ...) <- value

\S4method{colData}{FlexExperiment}(x, ...)

\S4method{colData}{FlexExperiment,ANY}(x, ...) <- value
}
\arguments{
\item{x, object}{A \code{\link{FlexExperiment}} object.}

\item{value}{An object of a class specified in the S4 method signature.}

\item{...}{Arguments passed to other metheds.}

\item{use.names}{Whether or not to set the dimension names of \code{x} as the row
names of output \code{\link[S4Vectors]{DataFrame}}.}
}
\value{
\itemize{
\item \code{length}: A single integer specifying the number of assays.
}

\itemize{
\item \code{names}: A character vector specifying the assay names.
}

\itemize{
\item \verb{names<-}: An updated \code{x} with new assay names.
}

\itemize{
\item \code{dim}: A two-element vector specifying the row and column numbers of
\code{x}.
}

\itemize{
\item \code{dimnames}: A list containing the row and column names of \code{x}.
}

\itemize{
\item \verb{rownames<-} and \verb{colnames<-}: An updated \code{x} with new dimension names.
}

\itemize{
\item \code{rowData}: A \code{\link[S4Vectors]{DataFrame}} specifying the meta data for
rows of \code{x}.
}

\itemize{
\item \verb{rowData<-}: An updated \code{x} with new metadata of rows.
}

\itemize{
\item \code{colData}: A \code{\link[S4Vectors]{DataFrame}} specifying the meta data for
columns of \code{x}.
}

\itemize{
\item \verb{colData<-}: An updated \code{x} with new metadata of columns.
}
}
\description{
A set of pre-defined accessor methods to extract information from a
\code{\link{FlexExperiment}}.
}

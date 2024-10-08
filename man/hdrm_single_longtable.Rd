% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdrm_single_methods.R
\name{hdrm_single_longtable}
\alias{hdrm_single_longtable}
\title{Test for one group high dimensional repeated measures}
\usage{
hdrm_single_longtable(
  data,
  hypothesis = "flat",
  value,
  subject,
  dimension,
  ...
)
}
\arguments{
\item{data}{a data.frame in longtable format.}

\item{hypothesis}{either "equal" or a quadratic numeric matrix.}

\item{value}{name or index of value column.}

\item{subject}{name or index of subject column.}

\item{dimension}{name or index of dimension column.}

\item{...}{further arguments are currently ignored.}
}
\value{
Returns a list with class "hdrm_single" with the components

\item{data}{the initial input data.}

\item{f}{the degrees of freedom \eqn{f}.}

\item{statistic}{the test statistic \eqn{W}.}

\item{tau}{the convergence parameter \eqn{\tau}.}

\item{H}{the hypothesis-matrix used.}

\item{hypothesis}{a character. Will be 'custom' if \code{hypothesis} was
given as a matrix, otherwise \code{hypothesis[1]}.}

\item{p.value}{the p-value of the test statistic.}

\item{dim}{a named vector giving the dimensions \eqn{d \times N} of \code{data}.}

\item{removed.cases}{number of subjects removed for having missing values.}
}
\description{
This function implements the methods outlined in
\insertCite{Pauly2015;textual}{hdrm} for data in a longtable format. For
data in a widetable format see \link{hdrm_single_widetable}.
}
\details{
The function can deal with missing values only in the \code{value} column.
Affected subjects are dropped with a warning. Missing values in any other
column of data will result in an error.

The test outlined in \insertCite{Pauly2015;textual}{hdrm} is performed for
the hypothesis given by \code{hypothesis}. The \code{hypothesis} can either be given as
a \code{character} or a \code{list}. The only legal character is "flat" (default).
"flat" stands for \eqn{H_0}: time profile is flat. Other
characters will result in an error.

Alternatively \code{hypothesis} can be a quadratic hypothesis matrix \eqn{T} with
the number of rows equal to the number of levels in the  of \code{data}. \eqn{T} must be
a projection matrix, meaning symmetrical and \eqn{T^2 = T}. A matrix that does not
match those criteria will result in an error.
}
\examples{
data(EEG)
names(EEG)
# name (index) of columns of interest are:
# value = value (2)
# subject = subject (4)
# dimension = dimension (7)

# call test
hdrm_single_longtable(data = EEG,
                      # test if time profiles are flat
                      hypothesis = "flat",
                      value = "value",
                      # columns can also be specified by their indices
                      subject = 4,
                      dimension = "dimension"
)


# hypothesis as list
hypothesis_matrix <- diag(40) - matrix(1/40, 40, 40)

hdrm_single_longtable(data = EEG,
                      # test if time profiles are flat
                      hypothesis = hypothesis_matrix,
                      value = "value",
                      subject = "subject",
                      dimension = "dimension"
)
}
\references{
\insertAllCited
}

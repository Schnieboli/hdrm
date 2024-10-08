% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdrm_grouped_methods.R
\name{hdrm_grouped_longtable}
\alias{hdrm_grouped_longtable}
\title{Test for multiple group high dimensional repeated measures}
\usage{
hdrm_grouped_longtable(
  data,
  hypothesis = c("whole", "sub", "interaction", "identical", "flat"),
  group,
  value,
  subject,
  dimension,
  subsampling = FALSE,
  B = "1000*N",
  ...
)
}
\arguments{
\item{data}{a data.frame in longtable format.}

\item{hypothesis}{either one of "whole", "sub" "interaction", "identical" or
"flat" or a named list with quadratic matrices \code{TW} and \code{TS} (see details).}

\item{group}{name or index of group column.}

\item{value}{name or index of value column.}

\item{subject}{name or index of subject column.}

\item{dimension}{name or index of dimension column.}

\item{subsampling}{\code{logical}. Specifying whether the subsampling versions for
all trace estimators should be used (see details).}

\item{B}{a \code{string} specifying a function of the number of subjects \eqn{N}.
Determines the number of subsamples used by the  subsampling trace
estimators.}

\item{...}{further arguments. Currently ignored.}
}
\value{
a named list of class "hdrm_grouped" with the components

\item{data}{the input data used.}

\item{f}{the degrees of freedom \eqn{f}.}

\item{tau}{the convergence parameter \eqn{\tau}.}

\item{H}{a named list with components \code{TW} and \code{TS} that give the
components of the hypothesis matrix.}

\item{hypothesis}{a character. Will be "custom" if \code{hypothesis} is a
list, otherwise \code{hypothesis[1]}.}

\item{p.value}{the \eqn{p}-value of the test statistic.}

\item{dim}{a named list with with number of factor levels \eqn{d} and
number of subjects \eqn{N} of \code{data}.}

\item{groups}{a named list with components number of groups \code{a} and
distribution of groups \code{table}.}

\item{removed.cases}{number of incomplete subjects removed.}

\item{subsamples}{number of subsamples used for subsampling
estimators.}

\insertNoCite{Sattler2018}{hdrm}
}
\description{
This function implements the methods outlined in
\insertCite{Sattler2018;textual}{hdrm} for data in longtable format. For
data in widetable format see \link{hdrm_grouped_widetable}.
}
\details{
The function can deal with missing values only in \code{value}. The test
is then performed without the affected subjects. Missing values in any of the
other columns of \code{data} will result in an error.

The test outlined in \insertCite{Sattler2018;textual}{hdrm} is performed for
the hypothesis \eqn{(\bm T_W \otimes \bm T_S)\bm\mu = \bm 0}, with \eqn{\bm
T_W}, \eqn{\bm T_S} given by \code{hypothesis}. The \code{hypothesis} can either be
given as a \code{character} or a \code{list}. Legal characters are "whole" (default),
"sub", "interaction", "identical" and "flat". The first three are given by (as
outlined in \insertCite{Sattler2018;textual}{hdrm}):
\itemize{
\item \eqn{H_0^{\overline W}\text{: } \left(\bm P_a \otimes \frac{1}{d}\bm J_d \right)\bm\mu=\bm 0}
(no main effect in whole plot factor group)
\item  \eqn{H_0^{\overline S}\text{: } \left(\frac{1}{a}\bm J_a \otimes \bm P_d \right)\bm\mu=\bm 0}
(no main effect in sub plot factor dimension)
\item \eqn{H_0^{WS}\text{: } \left(\bm P_a \otimes \bm P_d \right)\bm\mu=\bm 0}
(no interaction effect of whole plot and sub plot factor).
}
Additionally, "identical" and "flat" specify the hypotheses
\itemize{
\item  \eqn{H_0^{W}\text{: } \left(\bm P_a \otimes \bm I_d \right)\bm\mu=\bm 0}
(vectors of expexcted values are equal between all groups)
\item \eqn{H_0^{S}\text{: } \left(\bm I_a \otimes \bm P_d \right)\bm\mu=\bm 0}
(flat profile in all groups).
}
Other characters will result in an error.

Alternatively, \code{hypothesis} can be a named list with quadratic
matrices \code{TW} for the wholeplot-part and \code{TS} for the subplot-part of \eqn{\bm T
= \bm T_W \otimes \bm T_S}. If \code{hypothesis} is a list, \code{TW} and \code{TS} must be
projection matrices, meaning symmetrical and \eqn{\bm T_W^2 = \bm T_W},
\eqn{\bm T_S^2 = \bm T_S}. Differences up to the tolerance in \code{all.equal} are
ignored, while bigger differences will throw a warning. Also, the number of
rows and columns of \code{TW} must be equal to the number of groups and the number
of rows and columns of \code{TS} must be equal to the number of factor levels.
Lists that do not match those criteria will result in an error.

Note that for \code{subsampling = FALSE} results are still seed dependent,
because the computational heaviest trace estimator is still calculated using
subsamples. Also, depending on the choice of \code{B}, the non subsampling versions
might not be faster for small data.

That also means, that even for \code{subsampling = FALSE} the \code{p.value} depends
heavily on the choice of \code{B}, as it depends on the degrees of freedom \code{f},
which are always estimated by the subsampling version. For \code{subsampling = TRUE} the test statistic \eqn{W} is also seed dependent.

The number of subsamples \code{B} can also be a numeric value. However, this is not
advised, as \code{B} should be a function of \eqn{N} that goes to \eqn{\infty} for
\eqn{N\to\infty}.
}
\examples{
data(EEG)
names(EEG)
# name (index) of columns of interest are:
# group = group (1)
# value = value (2)
# subject = subject (4)
# dimension = dimension (7)

# call test
hdrm_grouped_longtable(data = EEG,
                       # test for no group effect
                       hypothesis = "whole",
                       group = "group",
                       value = "value",
                       # columns can also be specified by their indices
                       subject = 4,
                       dimension = "dimension",
                       subsampling = FALSE,
                       B = "100*N"
)


# test using all subsampling version of the trace estimators
hdrm_grouped_longtable(data = EEG,
                       # test for no time effect
                       hypothesis = "sub",
                       group = "group",
                       value = "value",
                       subject = "subject",
                       dimension = "dimension",
                       subsampling = TRUE,
                       # B can also be given as a number
                       B = 10000
)


# hypothesis as list
hypothesis_list <- list(TW = matrix(1/4, 4, 4),
                        TS = diag(40) - matrix(1/40, 40, 40)
) # equivalent to hypothesis = "sub"

# call test
hdrm_grouped_longtable(data = EEG,
                       hypothesis = hypothesis_list,
                       group = "group",
                       value = "value",
                       subject = "subject",
                       dimension = "dimension",
                       subsampling = FALSE,
                       B = "100*N"
)
rm(hypothesis_list)
}
\references{
\insertAllCited
}

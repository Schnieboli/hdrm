#' EEG measurements from 160 subjects
#'
#' A data frame containing four quantitative EEG variables measured at ten scalp
#' regions for each of 160 subjects. The resulting 40 repeated-measurement
#' dimensions are stored in long format. The data originate from a study
#' conducted by Höller et al. (2017) and are a part of the package `HRM`
#' described by Happ et al. (2018)
#'
#' @format A data frame with 6,400 rows and 7 variables:
#' \describe{
#'   \item{group}{Diagnostic group with levels `"SCC+"`, `"SCC-"`, `"MCI"`,
#'   and `"AD"`.}
#'   \item{value}{Numeric EEG-derived measurement.}
#'   \item{sex}{Recorded sex, encoded as `"M"` or `"W"`.}
#'   \item{subject}{Subject identifier.}
#'   \item{variable}{EEG variable coded from 1 to 4: activity, complexity,
#'   mobility, and brain rate.}
#'   \item{region}{Scalp region coded from 1 to 10: frontal, central,
#'   temporal, occipital, and parietal, each measured on the left and right.}
#'   \item{dimension}{Combined variable-region dimension coded from 1 to 40.}
#' }
#' 
#' The documentation was taken from Happ et al. (2018).
#'
#' @references Happ M, Harrar SW, Bathke AC (2018). “HRM: An R Package for Analysing High-dimensional Multi-factor Repeated Measures.” The R Journal, 10(1), 534–548. doi: 10.32614/RJ-2018-032.
#' @references Höller Y, Bathke AC, Uhl A, Strobl N, Lang A, Bergmann J, Nardone R, Rossini F, Zauner H, Kirschner M, Jahanbekam A, Trinka E, Staffen W (2017). “Combining SPECT and Quantitative EEG Analysis for the Automated Differential Diagnosis of Disorders with Amnestic Symptoms.” Frontiers in Aging Neuroscience, 9, 290. doi: 10.3389/fnagi.2017.00290.
#' @docType data
#' @keywords datasets
#' @name EEG
#' @usage EEG
"EEG"

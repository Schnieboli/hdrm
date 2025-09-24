#' EEG data of 160 subjects
#'
#' A dataset containing EEG data (Höller et al. 2017) of 160 subjects, 4 variables are measured at ten different locations.
#' The data set is part of the package HRM (Happ et al. 2017), which can be found at
#'  \url{https://cran.r-project.org/src/contrib/Archive/HRM/}.
#'
#' \itemize{
#'   \item `group`      Diagnostic group of the subject: Alzheimer's Disease (AD), Mild Cognitive Impairment (MCI), Subject Cognitive Complaints (SCC+, SCC-).
#'   \item `value`      Measured data of a subject at a specific variable and region.
#'   \item `sex`        Sex of the subject: Male (M) or Female (W).
#'   \item `subject`    A unique identification of a subject.
#'   \item `variable`   The variables measured are activity, complexity, mobility and brain rate coded from 1 to 4.
#'   \item `region`     Frontal left/right, central left/right, temporal left/right, occipital left/right, parietal left/right coded as 1 to 10.
#'   \item `dimension`  Mixing variable and region together, levels range from 1 to 40.
#' }
#' Documentation taken from Happ et al. (2017).
#' @references {
#'  Happ M, Harrar SW, Bathke AC (2017). “High-dimensional Repeated Measures.” Journal of Statistical Theory and Practice, 11(3), 468-477. doi:10.1080/15598608.2017.1307792.
#' }
#' @references {
#' Höller Y, Bathke AC, Uhl A, Strobl N, Lang A, Bergmann J, Nardone R, Rossini F, Zauner H, Kirschner M, Jahanbekam A, Trinka E, Staffen W (2017). “Combining SPECT and quantitative EEG analysis for the automated differential diagnosis of disorders with amnestic symptoms.” Front. Aging Neurosci., 9, 290. doi:10.3389/fnagi.2017.00290.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name EEG
#' @usage EEG
#' @format A data frame with 6400 rows and 7 columns.

"EEG"

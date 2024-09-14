#' EEG data of 160 subjects
#'
#' A dataset containing EEG data \insertCite{EEG_dataset}{hdrm} of 160 subjects, 4 variables are measured at ten different locations.
#'
#' \itemize{
#'   \item `group` Diagnostic group of the subject: Alzheimer's Disease (AD), Mild Cognitive Impairment (MCI), Subject Cognitive Complaints (SCC+, SCC-).
#'   \item `value` Measured data of a subject at a specific variable and region.
#'   \item `sex`` Sex of the subject: Male (M) or Female (W).
#'   \item `subject` A unique identification of a subject.
#'   \item `variable` The variables measured are activity, complexity, mobility and brain rate coded from 1 to 4.
#'   \item `region` Frontal left/right, central left/right, temporal left/right, occipital left/right, parietal left/right coded as 1 to 10.
#'   \item `dimension` Mixing variable and region together, levels range from 1 to 40.
#' }
#' Documentation taken from \insertCite{HRM_RJournal;textual}{hdrm}.
#' @references {
#'  \insertAllCited{}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name EEG
#' @usage EEG
#' @format A data frame with 6400 rows and 7 variables.

"EEG"

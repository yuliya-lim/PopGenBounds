#' Transform a list of matrices containing differentiation statistics into a tibble.
#' Exclude rows containing nan values.
#'
#' @param Diff_loci A list of matrices containing M, FST, G'ST and D values
#' for each loci in a sample. The values expected to be in the second column.
#'
#' @returns A tibble object without nan values.
#' @export
#'
#' @examples
#' library(tibble)
#' # Create example list of 2 matrices with 4 rows (M, FST, G'ST, D)
#' Diff_loci <- list(
#'   matrix(c("M", 0.85, NA, NA,
#'            "FST", 0.01, NA, NA,
#'            "G'ST", 0.02, NA, NA,
#'            "D", 0.015, NA, NA),
#'          ncol = 4, byrow = TRUE,
#'          dimnames = list(NULL, c("statistic", "value", "bound", "bound_closeness"))),
#'
#'   matrix(c("M", 0.5, NA, NA,
#'            "FST", 0.02, NA, NA,
#'            "G'ST", 0.03, NA, NA,
#'            "D", 0.025, NA, NA),
#'          ncol = 4, byrow = TRUE,
#'          dimnames = list(NULL, c("statistic", "value", "bound", "bound_closeness")))
#' )
#'
#' # Convert to tibble and remove rows with NaNs
#' clean_tibble <- turn_diff_in_tibble(Diff_loci)
#' print(clean_tibble)
turn_diff_in_tibble <- function(Diff_loci){
  diff.tib = tibble(M=unlist(sapply(Diff_loci,function(x){x[1,2]})),
                     FST=unlist(sapply(Diff_loci,function(x){x[2,2]})),
                     GpST=unlist(sapply(Diff_loci,function(x){x[3,2]})),
                     D=unlist(sapply(Diff_loci,function(x){x[4,2]})))

  diff.tib_clean <- diff.tib[!apply(is.nan(as.matrix(diff.tib)), 1, any), ]

  return(diff.tib_clean)
}

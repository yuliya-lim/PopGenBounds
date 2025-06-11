#' Compute differentiation statistics in case of K=3 subpopulations
#'
#' @param list_loci A list of K × 2 matrices of allele frequencies, one per locus.
#' Each matrix contains frequencies of the reference and alternative alleles (columns)
#' for each subpopulation (rows).
#'
#' @returns Four lists of matrices containing differentiation statistics FST, G'ST, D and a vector of frequencies M.
#' The first list corresponds to 3 subpopulations together, the remaining lists
#' correspond to pair-wise differentiation statistics for each pair of subpopulations out of 3 subpopulations.
#' @export
#'
#' @examples
#' # Toy example with 2 loci and 3 subpopulations
#' locus1 <- matrix(c(0.8, 0.2,   # Subpop 1
#'                    0.6, 0.4,   # Subpop 2
#'                    0.7, 0.3),  # Subpop 3
#'                  nrow = 3, byrow = TRUE)
#'
#' locus2 <- matrix(c(0.5, 0.5,
#'                    0.4, 0.6,
#'                    0.6, 0.4),
#'                  nrow = 3, byrow = TRUE)
#'
#' list_loci <- list(locus1, locus2)
#'
#' # Compute differentiation statistics
#' diff_stats <- compute_Diff_3subpop(list_loci)
#'
#' # Access 3-subpop differentiation stats for first locus
#' print(diff_stats$D_123[[1]])
#'
#' # Access pairwise differentiation (subpop 1 vs 2) for second locus
#' print(diff_stats$D_12[[2]])
compute_Diff_3subpop <- function(list_loci){
  # check whether K=3 (to do)

  # Extract couple of rows (subpops) from 3-subpop matrices
  list_loci_12 <- lapply(list_loci, function(mat) mat[c(1, 2), ])
  list_loci_23 <- lapply(list_loci, function(mat) mat[c(2, 3), ])
  list_loci_13 <- lapply(list_loci, function(mat) mat[c(1, 3), ])

  # Compute stats (F_ST, G_ST, D)
  D_loci_123 = lapply(list_loci, Diff)
  D_loci_12 = lapply(list_loci_12, Diff)
  D_loci_23 = lapply(list_loci_23, Diff)
  D_loci_13 = lapply(list_loci_13, Diff)

  return(list(D_123=D_loci_123, D_12=D_loci_12, D_13=D_loci_13, D_23=D_loci_23))

}

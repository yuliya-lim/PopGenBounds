#' Produces a list of allele frequency matrices for each locus,
#' formatted for input to the Diff() function
#'
#' @param data_frequencies A data frame of allele frequencies with dimensions (L × K),
#' where L is the number of loci (rows), and K is the number of subpopulations (columns)
#'
#' @returns A list of K × 2 matrices of allele frequencies, one per locus.
#' Each matrix contains frequencies of the reference and alternative alleles
#' for each subpopulation.
#' @export
#'
#' @examples
#' # Create a toy dataframe of allele frequencies for 3 loci and 2 subpopulations
#' df <- data.frame(
#'   Subpop1 = c(0.8, 0.6, 0.4),
#'   Subpop2 = c(0.7, 0.5, 0.3)
#' )
#' # Use the function to generate a list of 2x2 frequency matrices
#' freq_list <- make_popgen_input(df)
make_popgen_input <- function(data_frequencies) {
  list_loci <- lapply(1:nrow(data_frequencies), function(i){
    mat <- as.matrix(data_frequencies[i,])
    second_allele <- 1 - mat[1, ]  # Compute 1-x for each element in the first row
    new_mat <- rbind(mat, second_allele)
    rownames(new_mat) <- c("First", "Second")
    return(t(new_mat))
  }
  )
  return(list_loci)
}

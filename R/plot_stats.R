#' Produces a combined plot with raw and normalized differentiation statistics
#'
#' @param Diff_loci A list of matrices containing values of differentiation statistics FST, G'ST, D
#' and a vector of frequencies M.
#' @param K Number of subpopulations
#' @param title Title of the produced plot
#'
#' @returns A ggplot object
#' @export
#'
#' @examples
#' # Example with 2 loci and 2 subpopulations
#' library(tibble)
#' library(patchwork)
#' library(ggplot2)
#' library(glue)
#' locus1 <- matrix(c(0.8, 0.2,    # Subpop 1
#'                    0.6, 0.4),   # Subpop 2
#'                  nrow = 2, byrow = TRUE)
#'
#' locus2 <- matrix(c(0.5, 0.5,
#'                    0.4, 0.6),
#'                  nrow = 2, byrow = TRUE)
#'
#' list_loci <- list(locus1, locus2)
#' Diff_loci <- lapply(list_loci, Diff)
#' combined_plot <- plot_stats(Diff_loci, K=2, title=glue::glue("K=2: Subpop1, Subpop2"))
plot_stats <- function(Diff_loci, K=2, title="") {

  lcnec.tib_clean <- turn_diff_in_tibble(Diff_loci)

  gglcnec = ggbounds_new(M=lcnec.tib_clean$M,FST=lcnec.tib_clean$FST,GpST=lcnec.tib_clean$GpST,D=lcnec.tib_clean$D,K=K)

  combined_plot <-
    (gglcnec[[1]][[1]] + gglcnec[[2]][[1]]) /
    (gglcnec[[1]][[2]] + gglcnec[[2]][[2]]) /
    (gglcnec[[1]][[3]] + gglcnec[[2]][[3]]) +
    patchwork::plot_annotation(title = title,
                    theme = theme(
                      plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 20))
                    )
    )

  return(combined_plot)
}

#' Plot differentiation statistics in case of K = 3 subpopulations
#'
#' @param D_loci Four lists of matrices containing differentiation statistics FST, G'ST, D and a vector of frequencies M.
#' The first list corresponds to 3 subpopulations together, the remaining lists
#' correspond to pair-wise differentiation statistics for each pair of subpopulations out of 3 subpopulations.
#' @param subpop_names A string vector containing names of subpopulations to be plotted.
#'
#' @returns A ggplot object.
#' @import glue
#' @export
#'
#' @examples
#' library(patchwork)
#' library(glue)
#' library(tibble)
#' # Initialise allele frequencies for 3 subpopulations in 2 loci
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
#' # Compute differentiation statistics
#' diff_stats <- compute_Diff_3subpop(list(locus1, locus2))
#'
#' # Produce plot of differentiation statistics values
#' final_plot <- plot_stats_3subpop(diff_stats, c("Subpop1", "Subpop2", "Subpop3"))
plot_stats_3subpop <- function(D_loci, subpop_names){
  plot123 <- plot_stats(D_loci$D_123, K=3, title=glue("K=3: {subpop_names[1]}, {subpop_names[2]}, {subpop_names[3]}"))
  plot12 <- plot_stats(D_loci$D_12, K=2, title=glue("K=2: {subpop_names[1]}, {subpop_names[2]}"))
  plot23 <- plot_stats(D_loci$D_23, K=2, title=glue("K=2: {subpop_names[2]}, {subpop_names[3]}"))
  plot13 <- plot_stats(D_loci$D_13, K=2, title=glue("K=2: {subpop_names[1]}, {subpop_names[3]}"))

  final_plot <- wrap_elements(plot123) /
    wrap_elements(plot12) /
    wrap_elements(plot23) /
    wrap_elements(plot13)

  return(final_plot)
}

#' Pipeline to compute distributions of differnetiation statistics for LNEN data
#'
#' @param data_frame_list A list of dataframes with LNEN samples
#' @param K Number of subpopulations. K should be the same for all samples in data_frame_list
#'
#' @export
run_lnen_plotting <- function(data_frame_list, K=2) {
  for (df in data_frame_list) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
    data_clean <- filter_data(df, K)
    list_freq <- make_popgen_input(data_clean[,3:(2+K)])
    if (K==2){
      Diff_loci <- lapply(list_freq, Diff)
      combined_plot <- plot_stats(Diff_loci, K, title=glue("K=2: {subpop_names[1]}, {subpop_names[2]}"))
    }
    else if (K==3){
      D_loci = compute_Diff_3subpop(list_freq)
      combined_plot <- plot_stats_3subpop(D_loci, subpop_names)
      combined_plot
    }
    save_plots(subpop_names, combined_plot, type="stats", K)
  }
}

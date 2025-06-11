#' Computes mean values and normalized mean valuesof FST, G'ST and D
#' statistics over all loci in a given sample.
#'
#' @param Diff_loci A list of matrices containing values of differentiation statistics FST, G'ST, D
#' and a vector of frequencies M.
#' @param K Number of subpopulations
#'
#' @return A named list with six numeric values:
#' \describe{
#'   \item{FST_mean}{Mean raw FST value across loci}
#'   \item{GST_mean}{Mean raw G'ST value across loci}
#'   \item{D_mean}{Mean raw Jost's D value across loci}
#'   \item{FST_norm_mean}{Mean normalized FST value across loci}
#'   \item{GST_norm_mean}{Mean normalized G'ST value across loci}
#'   \item{D_norm_mean}{Mean normalized Jost's D value across loci}
#' }
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
#' mean_stats <- compute_mean_stats(Diff_loci, K=2)
compute_mean_stats <- function(Diff_loci, K){

  lcnec.tib_clean <- turn_diff_in_tibble(Diff_loci)

  FST_mean <- mean(lcnec.tib_clean$FST, na.rm = TRUE)
  FST_norm <- lcnec.tib_clean$FST / sapply(lcnec.tib_clean$M, function(m) Fup(K, m))
  FST_norm_mean <- mean(FST_norm, na.rm=T)

  GST_mean <- mean(lcnec.tib_clean$GpST, na.rm = TRUE)
  GST_norm <- lcnec.tib_clean$GpST / sapply(lcnec.tib_clean$M, function(m) Gpup(K, m))
  GST_norm_mean <- mean(GST_norm, na.rm=T)

  D_mean <- mean(lcnec.tib_clean$D, na.rm = TRUE)
  D_norm <- lcnec.tib_clean$D / sapply(lcnec.tib_clean$M, function(m) Dup(K, m))
  D_norm_mean <- mean(D_norm, na.rm=T)

  return(list(
    FST_mean = FST_mean,
    GST_mean = GST_mean,
    D_mean = D_mean,
    FST_norm_mean = FST_norm_mean,
    GST_norm_mean = GST_norm_mean,
    D_norm_mean = D_norm_mean
  ))
}

#' Computes mean and normalized mean values of FST, G'ST, D for a list of samples.
#'
#' @param data_frame_list A list of dataframes containing allele fraquencies for loci (rows)
#' in multiple subpopulations (columns)
#' @param K Number of subpopulations
#'
#' @returns A list of named lists containing six numeric values:
#' \describe{
#'   \item{FST_mean}{Mean raw FST value across loci}
#'   \item{GST_mean}{Mean raw G'ST value across loci}
#'   \item{D_mean}{Mean raw Jost's D value across loci}
#'   \item{FST_norm_mean}{Mean normalized FST value across loci}
#'   \item{GST_norm_mean}{Mean normalized G'ST value across loci}
#'   \item{D_norm_mean}{Mean normalized Jost's D value across loci}
#' } for a provided list of samples
#'
#' @export
get_mean_stats <- function(data_frame_list, K=2){
  mean_stats_list <- list()
  for (df in data_frame_list) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
    sample <- subpop_names[[1]]
    data_clean <- filter_data(df, K)
    list_freq <- make_popgen_input(data_clean[,3:(2+K)])
    if (K==2){
      Diff_loci <- lapply(list_freq, Diff)
      mean_stats <- compute_mean_stats(Diff_loci, K)
    }
    else if (K==3){
      print("Computing Diff for all pairs of subpopulations")
      D_loci = compute_Diff_3subpop(list_freq)
      mean_stats <- compute_mean_stats(D_loci$D_123, K)
    }
    mean_stats_list[[sample]] <- mean_stats
  }
  return(mean_stats_list)
}

#' Plots mean values of FST, G'ST and D and their normalised mean values for a given list of samples.
#'
#' @param mean_stats_list A list of named list containing six numeric values:
#' \describe{
#'   \item{FST_mean}{Mean raw FST value across loci}
#'   \item{GST_mean}{Mean raw G'ST value across loci}
#'   \item{D_mean}{Mean raw Jost's D value across loci}
#'   \item{FST_norm_mean}{Mean normalized FST value across loci}
#'   \item{GST_norm_mean}{Mean normalized G'ST value across loci}
#'   \item{D_norm_mean}{Mean normalized Jost's D value across loci}
#' }
#' Each named list corresponds to a different data sample.
#'
#' @returns A ggplot object
#' @export
#'
plot_mean_stats <- function(mean_stats_list){
  mean_stats_df <- bind_rows(lapply(mean_stats_list, as_tibble), .id = "Sample")
  sample_order <- c("SINET8M", "SINET9M", "LNET6T", "LNET10T", "LCNEC3T", "LCNEC4T", "PANEC1T")
  mean_stats_df$Sample <- factor(mean_stats_df$Sample, levels = sample_order)

  # Step 1: Reshape to long format
  mean_stats_long <- mean_stats_df %>%
    select(Sample, FST_mean, GST_mean, D_mean) %>%
    pivot_longer(cols = -Sample, names_to = "Metric", values_to = "Value")

  mean_norm_stats_long <- mean_stats_df %>%
    select(Sample, FST_norm_mean, GST_norm_mean, D_norm_mean) %>%
    pivot_longer(cols = -Sample, names_to = "Metric", values_to = "Value")

  # Define the order of metrics as you want them in the legend
  metric_order <- c("FST_mean", "GST_mean", "D_mean")
  metric_norm_order <- c("FST_norm_mean", "GST_norm_mean", "D_norm_mean")

  # Set factor levels to control legend order
  mean_stats_long$Metric <- factor(mean_stats_long$Metric, levels = metric_order)
  mean_norm_stats_long$Metric <- factor(mean_norm_stats_long$Metric, levels = metric_norm_order)

  # Step 2: Plot
  g1 <- ggplot(mean_stats_long, aes(x = Sample, y = Value, color = Metric, group = Metric)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    theme_bw() +
    labs(title = "Mean Diff Statistics per Sample",
         x = "Sample", y = "Mean Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  g2 <- ggplot(mean_norm_stats_long, aes(x = Sample, y = Value, color = Metric, group = Metric)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    theme_bw() +
    labs(title = "Mean Normalised Diff Statistics per Sample",
         x = "Sample", y = "Mean Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  combined_plot <- g1 + g2
  return(combined_plot)
}

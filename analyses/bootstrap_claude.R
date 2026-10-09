# ---------------------------------------------------------------
# bootstrap.R
# Bootstrap confidence intervals for statistics computed from
# allele frequency matrices (Fst, G'st, D, ...)
#
# Resampling unit: individuals, within each population, with replacement.
#
#   res <- boot_stat(toad_genotypes, stat_fun = my_fst, B = 1000)
#   res$estimate; res$ci
# ---------------------------------------------------------------


#' Bootstrap a statistic over individuals within populations
#'
#' Resamples individuals with replacement within each population (keeping
#' the original sample sizes), recomputes the allele frequency matrices and
#' the statistic for each replicate, and returns percentile confidence
#' intervals. The allele set of the original data is used for every
#' replicate, so alleles not drawn get frequency 0.
#'
#' @param genotypes Table from [microsat_genotypes()].
#' @param stat_fun Function taking a list of allele frequency matrices (as
#'   returned by [microsat_freq()]) and returning a number, a vector or a
#'   matrix (e.g. pairwise values).
#' @param B Number of bootstrap replicates.
#' @param level Confidence level of the intervals.
#' @param ci_cols For a statistic returning a data frame/tibble: columns to
#'   compute intervals for (default: all numeric columns).
#' @param ... Further arguments passed to `stat_fun`.
#'
#' @return A list with
#'   \describe{
#'     \item{}{If `stat_fun` returns a data frame/tibble: `estimate` is that
#'       table with added columns `<col>_lower` and `<col>_upper`, and `boot`
#'       stacks all replicates with a `replicate` column. Otherwise:}
#'     \item{estimate}{Statistic computed on the original data.}
#'     \item{ci}{Percentile confidence interval(s), see [boot_ci()].}
#'     \item{boot}{Bootstrap values: a vector (scalar statistic) or an array
#'       whose last dimension indexes the replicates.}
#'     \item{B, level}{Settings used.}
#'   }
#' @export
#'
#' @examples
#' # mean expected heterozygosity, as a simple example statistic
#' mean_he <- function(fl) mean(sapply(fl, function(m) mean(1 - rowSums(m^2), na.rm = TRUE)))
#' set.seed(1)
#' res <- boot_stat(toad_genotypes, mean_he, B = 50)
#' res$estimate
#' res$ci
boot_stat <- function(genotypes, stat_fun, B = 1000, level = 0.95, ci_cols = NULL, ...) {
  lv <- microsat_allele_levels(genotypes)              # fixed allele set (original data)
  estimate <- stat_fun(microsat_freq(genotypes, lv), ...)   # stats of original data

  boot <- lapply(seq_len(B), function(b) {            # bootstrap distribution of stats
    idx <- boot_indices(genotypes$pop)
    stat_fun(microsat_freq(genotypes[idx, ], lv), ...)
  })
  if (is.data.frame(estimate)) {                       # tibble output: CI per numeric column
    return(list(estimate = boot_ci_df(estimate, boot, level, ci_cols),
                boot = do.call(rbind, Map(cbind, replicate = seq_len(B), boot)),
                B = B, level = level))
  }

  boot <- simplify2array(boot, higher = TRUE)          # vector, or array with replicates last
  list(estimate = estimate, ci = boot_ci(boot, level), boot = boot, B = B, level = level)
}


# Internal: add <col>_lower / <col>_upper to a data frame of estimates.
# Assumes every replicate returns the same rows in the same order.
boot_ci_df <- function(estimate, boot, level, ci_cols = NULL) {
  if (is.null(ci_cols)) ci_cols <- names(estimate)[vapply(estimate, is.numeric, logical(1))]
  arr  <- simplify2array(lapply(boot, function(d) as.matrix(d[ci_cols])), higher = TRUE)
  lims <- boot_ci(arr, level)                          # rows x ci_cols matrices
  estimate[paste0(ci_cols, "_lower")] <- as.data.frame(lims$lower)
  estimate[paste0(ci_cols, "_upper")] <- as.data.frame(lims$upper)
  estimate
}


#' Percentile confidence intervals from bootstrap values
#'
#' @param boot Bootstrap values: a vector, or an array whose last dimension
#'   indexes the replicates (as in the `boot` element of [boot_stat()]).
#' @param level Confidence level.
#'
#' @return For a vector, the lower and upper limits. For an array, a list
#'   with `lower` and `upper`, each with the shape of one replicate
#'   (e.g. a population x population matrix).
#' @export
boot_ci <- function(boot, level = 0.95) {
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  if (is.null(dim(boot))) return(stats::quantile(boot, probs, na.rm = TRUE))

  margins <- seq_len(length(dim(boot)) - 1)            # all dimensions except replicates
  list(lower = apply(boot, margins, stats::quantile, probs = probs[1], na.rm = TRUE),
       upper = apply(boot, margins, stats::quantile, probs = probs[2], na.rm = TRUE))
}


# Internal: row indices of one bootstrap sample, resampling individuals
# with replacement within each population (sample sizes unchanged).
boot_indices <- function(pop) {
  unlist(lapply(split(seq_along(pop), pop),
                function(ix) ix[sample.int(length(ix), replace = TRUE)]),
         use.names = FALSE)
}

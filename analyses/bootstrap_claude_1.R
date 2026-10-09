# ---------------------------------------------------------------
# bootstrap.R
# Bootstrap confidence intervals for per-locus statistics (Fst, G'st, D)
#
# Resampling unit: individuals, within each population, with replacement.
#
#   res <- boot_stat(toad_genotypes, stat_fun = Diff, B = 1000)
#   res$estimate[[1]]   # locus 1: statistic, value, bound, ..., value_lower, value_upper
# ---------------------------------------------------------------


#' Bootstrap confidence intervals for per-locus statistics
#'
#' Resamples individuals with replacement within each population (keeping
#' the original sample sizes), recomputes the allele frequency matrices and
#' applies `stat_fun` to the matrix of each locus. Percentile confidence
#' intervals are added to the estimates of the original data. The allele set
#' of the original data is used in every replicate, so alleles not drawn get
#' frequency 0.
#'
#' @param genotypes Table from [microsat_genotypes()].
#' @param stat_fun Function applied to one allele frequency matrix
#'   (populations x alleles), returning a tibble with one row per statistic
#'   (e.g. [Diff()]). Rows must come out in the same order every time.
#' @param B Number of bootstrap replicates.
#' @param level Confidence level of the intervals.
#' @param ci_cols Columns of the `stat_fun` output to compute intervals for.
#' @param keep_boot If `TRUE`, also return the bootstrap values.
#' @param ... Further arguments passed to `stat_fun`.
#'
#' @return A list with
#'   \describe{
#'     \item{estimate}{List (one element per locus) of the `stat_fun` output
#'       for the original data, with added columns `<col>_lower` and
#'       `<col>_upper` for each column in `ci_cols`.}
#'     \item{boot}{If `keep_boot = TRUE`: list (one element per locus) of the
#'       bootstrap outputs stacked into one table with a `replicate` column;
#'       otherwise `NULL`.}
#'     \item{B}{Number of bootstrap replicates used.}
#'     \item{level}{Confidence level of the intervals.
#'   }
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(2021)
#' res <- boot_stat(toad_genotypes, stat_fun = Diff, B = 20)
#' res$estimate[[1]]
#' }
boot_stat <- function(genotypes, stat_fun, B = 1000, level = 0.95,
                      ci_cols = "value", keep_boot = FALSE, ...) {
  lv <- microsat_allele_levels(genotypes)                 # fixed allele set (original data)
  freq_orig <- microsat_freq(genotypes, lv)
  estimate  <- lapply(freq_orig, stat_fun, ...)           # stats of original data, per locus

  boot <- lapply(seq_len(B), function(b) {                # bootstrap distribution of stats
    idx <- boot_indices(genotypes$pop)
    freq_boot <- microsat_freq(genotypes[idx, ], lv)
    lapply(freq_boot, stat_fun, ...)
  })

  # confidence intervals, locus by locus
  result <- lapply(seq_along(estimate), function(k) {
    boot_k <- lapply(boot, function(r) r[[k]])            # locus k from every replicate
    boot_ci_df(estimate[[k]], boot_k, level, ci_cols)
  })
  names(result) <- names(estimate)

  # optional: bootstrap values, stacked per locus
  boot_tables <- NULL
  if (keep_boot) {
    boot_tables <- lapply(seq_along(estimate), function(k)
      do.call(rbind, lapply(seq_len(B), function(b) cbind(replicate = b, boot[[b]][[k]]))))
    names(boot_tables) <- names(estimate)
  }

  list(estimate = result, boot = boot_tables, B = B, level = level)
}


# Internal: add <col>_lower / <col>_upper (percentile limits) to the
# estimate table of one locus. boot_k = list of B bootstrap tables of that
# locus; intervals are computed row by row (same row order in every table).
boot_ci_df <- function(estimate, boot_k, level, ci_cols) {
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  for (col in ci_cols) {
    # rows x B matrix: column b holds replicate b
    vals <- matrix(unlist(lapply(boot_k, function(d) d[[col]])), nrow = nrow(estimate))
    estimate[[paste0(col, "_lower")]] <- apply(vals, 1, stats::quantile, probs = probs[1], na.rm = TRUE)
    estimate[[paste0(col, "_upper")]] <- apply(vals, 1, stats::quantile, probs = probs[2], na.rm = TRUE)
  }
  estimate
}


# Internal: row indices of one bootstrap sample, resampling individuals
# with replacement within each population (sample sizes unchanged).
boot_indices <- function(pop) {
  unlist(lapply(split(seq_along(pop), pop),
                function(ix) ix[sample.int(length(ix), replace = TRUE)]),
         use.names = FALSE)
}

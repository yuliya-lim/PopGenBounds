# ---------------------------------------------------------------
# pairwise.R
# Pairwise differentiation statistics (K = 2) for all population pairs:
# per locus, and their arithmetic mean over loci, with bootstrap CIs.
# Output mirrors boot_stat(): res$estimate (one table per locus) and res$mean.
#
#   pw  <- pairwise_stats(microsat_freq(toad_genotypes))   # observed values
#   res <- boot_pairwise(toad_genotypes, B = 1000)          # with CIs
# ---------------------------------------------------------------


#' Pairwise differentiation statistics
#'
#' For every pair of populations, computes \eqn{M}, \eqn{F_{ST}},
#' \eqn{G'_{ST}} and \eqn{D} at each locus with the same formulas as [Diff()]
#' applied to the two rows of the pair (\eqn{K = 2}), and their arithmetic
#' mean over loci. Loci where a statistic is undefined for a pair (e.g. both
#' populations fixed for the same allele) are skipped in the mean.
#'
#' @param freq_list List of allele frequency matrices (populations x alleles),
#'   one per locus, as returned by [microsat_freq()].
#'
#' @return A list with
#'   \describe{
#'     \item{estimate}{List (one element per locus) of tables with columns
#'       `pop1`, `pop2`, `statistic` and `value`, one row per pair and statistic.}
#'     \item{mean}{Table with the same columns: mean over loci of each statistic.}
#'   }
#' @export
pairwise_stats <- function(freq_list) {
  pops  <- rownames(freq_list[[1]])
  up    <- which(upper.tri(diag(length(pops))), arr.ind = TRUE)     # array of index pairs: each pair of pop once
  pairs <- data.frame(pop1 = pops[up[, 1]], pop2 = pops[up[, 2]])
  v     <- pairwise_values(freq_list, up)   # returns 'per_locus' (3-dim array pairs x stats x loci) and 'mean' (pairs x stats array of means over loci)

  estimate <- lapply(seq_along(freq_list), function(l)
    ci_table(matrix(v$per_locus[, , l], nrow = nrow(up)), NULL, pairs))
  names(estimate) <- names(freq_list)
  list(estimate = estimate, mean = ci_table(v$mean, NULL, pairs))
}


#' Bootstrap confidence intervals for pairwise statistics
#'
#' Resamples individuals with replacement within each population (as in
#' [boot_stat()]); each replicate resamples all populations once, and every
#' pair is computed from it. Confidence intervals, bias and bias-corrected
#' values are computed for the per-locus statistics and for their mean over
#' loci.
#'
#' @param genotypes Table from [microsat_genotypes()].
#' @param B Number of bootstrap replicates.
#' @param level Confidence level.
#' @param type `"basic"` (default) or `"percentile"` intervals.
#' @param locus_ci If `TRUE` (default), compute intervals for each locus as
#'   well as for the mean. `FALSE` gives per-locus values without intervals
#'   and uses much less memory.
#'
#' @return A list with
#'   \describe{
#'     \item{estimate}{List (one element per locus) of tables with columns
#'       `pop1`, `pop2`, `statistic`, `value` (original data) and, if
#'       `locus_ci = TRUE`, `value_lower`, `value_upper`, `value_bias` and
#'       `value_corrected`.}
#'     \item{mean}{Same table for the mean over loci, always with intervals.}
#'     \item{B}{Number of bootstrap replicates used.}
#'     \item{level}{Confidence level of the intervals.}
#'     \item{type}{Interval type used.}
#'   }
#' @export
boot_pairwise <- function(genotypes, B = 1000, level = 0.95,
                          type = c("basic", "percentile"), locus_ci = TRUE) {
  type <- match.arg(type)
  genotypes$pop <- droplevels(genotypes$pop)
  lv    <- microsat_allele_levels(genotypes)
  freq  <- microsat_freq(genotypes, lv)
  pops  <- rownames(freq[[1]])
  up    <- which(upper.tri(diag(length(pops))), arr.ind = TRUE)     # each pair once
  pairs <- data.frame(pop1 = pops[up[, 1]], pop2 = pops[up[, 2]])
  P <- nrow(up)                                                     # number of pairs
  L <- length(lv)                                                   # number of loci

  obs <- pairwise_values(freq, up)                                  # original data

  # bootstrap values: pairs x statistics x (loci) x replicates
  boot_mean <- array(NA_real_, dim = c(P, 4, B))
  boot_loc  <- if (locus_ci) array(NA_real_, dim = c(P, 4, L, B)) else NULL
  for (b in seq_len(B)) {
    if (b %% 100 == 0) message("bootstrap ", b, "/", B)
    idx <- boot_indices(genotypes$pop)
    v   <- pairwise_values(microsat_freq(genotypes[idx, ], lv), up)
    boot_mean[, , b] <- v$mean
    if (locus_ci) boot_loc[, , , b] <- v$per_locus
  }

  estimate <- lapply(seq_len(L), function(l) {
    vals <- if (locus_ci) array(boot_loc[, , l, ], dim = c(P, 4, B)) else NULL
    ci_table(matrix(obs$per_locus[, , l], nrow = P), vals, pairs, level, type)
  })
  names(estimate) <- names(lv)

  list(estimate = estimate,
       mean     = ci_table(obs$mean, boot_mean, pairs, level, type),
       B = B, level = level, type = type)
}

# ---- internal helpers ---------------------------------------------


# Internal: HS, HT and M for all pairs of populations at one locus.
# f = populations x alleles frequency matrix. For pair (i, j):
#   HS = mean of the two within-population heterozygosities
#   HT = 1 - sum of squared mean frequencies of the pair
#   M  = largest mean frequency of the pair
pairwise_H <- function(f) {
  s  <- rowSums(f^2)                                          # sum of squared freqs per population
  # Compute HS from formula: Hs = 1 - S/K, with S = sum of squared frequencies over alleles and populations
  # (for each pair of populations: HS[i,j] = HS(pop_i,pop_j))).
  HS <- 1 - outer(s, s, "+") / 2
  #(for each pair of populations: HT[i,j] = HT(pop_i,pop_j))).
  HT <- 1 - (outer(s, s, "+") + 2 * tcrossprod(f)) / 4        # sum(((f_i + f_j) / 2)^2), all pairs at once.
  # Compute M - maximal mean allele frequency over populations.
  M  <- matrix(0, nrow(f), nrow(f))
  for (a in seq_len(ncol(f))) M <- pmax(M, outer(f[, a], f[, a], "+") / 2)
  list(HS = HS, HT = HT, M = M)
}

# M, FST, G'ST and D for all pairs at one locus (Diff() formulas with K = 2).
# Returns a list of four population x population matrices.
pairwise_locus <- function(f) {
  h   <- pairwise_H(f)
  FST <- ifelse(h$HT < 1e-12, 0, 1 - h$HS / h$HT)
  list(M      = h$M,
       FST    = FST,
       "G'ST" = FST * (1 + h$HS) / (1 - h$HS),
       D      = 2 * (h$HT - h$HS) / (1 - h$HS))
}

# Per-locus values and their mean over loci, for the pairs in `up`.
# Returns per_locus (pairs x 4 statistics x loci) and mean (pairs x 4 statistics).
pairwise_values <- function(freq_list, up) {
  per_locus <- lapply(freq_list, function(f) {
    s <- pairwise_locus(f)
    matrix(sapply(s, function(m) m[up]), nrow = nrow(up))     # pairs x statistics
  })
  arr  <- array(unlist(per_locus), dim = c(nrow(up), 4, length(freq_list)))   # 3-dimensional array, with loci as the third dim
  avg <- matrix(apply(arr, c(1, 2), mean, na.rm = TRUE), nrow = nrow(up))
  list(per_locus = arr, mean = avg)
}

# Long table (one row per pair and statistic) from estimates and, optionally,
# bootstrap values. est = pairs x 4 statistics; vals = pairs x 4 x B or NULL.
ci_table <- function(est, vals, pairs, level = 0.95, type = "basic") {
  stat_names <- c("M", "FST", "G'ST", "D")
  P <- nrow(pairs)
  out <- lapply(seq_along(stat_names), function(k) {
    tab <- data.frame(pairs, statistic = stat_names[k], value = est[, k])
    if (!is.null(vals)) {
      v     <- matrix(vals[, k, ], nrow = P)                    # pairs x B
      probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
      q_lo  <- apply(v, 1, stats::quantile, probs = probs[1], na.rm = TRUE)
      q_hi  <- apply(v, 1, stats::quantile, probs = probs[2], na.rm = TRUE)
      if (type == "percentile") {
        tab$value_lower <- q_lo
        tab$value_upper <- q_hi
      } else {                                                  # basic: reflected around the estimate
        tab$value_lower <- 2 * tab$value - q_hi
        tab$value_upper <- 2 * tab$value - q_lo
      }
      tab$value_bias      <- rowMeans(v, na.rm = TRUE) - tab$value
      tab$value_corrected <- tab$value - tab$value_bias
    }
    tab
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

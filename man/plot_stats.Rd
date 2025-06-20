% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_stats.R
\name{plot_stats}
\alias{plot_stats}
\title{Produces a combined plot with raw and normalized differentiation statistics}
\usage{
plot_stats(Diff_loci, K = 2, title = "")
}
\arguments{
\item{Diff_loci}{A list of matrices containing values of differentiation statistics FST, G'ST, D
and a vector of frequencies M.}

\item{K}{Number of subpopulations}

\item{title}{Title of the produced plot}
}
\value{
A ggplot object
}
\description{
Produces a combined plot with raw and normalized differentiation statistics
}
\examples{
# Example with 2 loci and 2 subpopulations
library(tibble)
library(patchwork)
library(ggplot2)
library(glue)
locus1 <- matrix(c(0.8, 0.2,    # Subpop 1
                   0.6, 0.4),   # Subpop 2
                 nrow = 2, byrow = TRUE)

locus2 <- matrix(c(0.5, 0.5,
                   0.4, 0.6),
                 nrow = 2, byrow = TRUE)

list_loci <- list(locus1, locus2)
Diff_loci <- lapply(list_loci, Diff)
combined_plot <- plot_stats(Diff_loci, K=2, title=glue::glue("K=2: Subpop1, Subpop2"))
}

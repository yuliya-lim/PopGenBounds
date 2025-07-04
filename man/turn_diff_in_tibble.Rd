% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{turn_diff_in_tibble}
\alias{turn_diff_in_tibble}
\title{Transform a list of matrices containing differentiation statistics into a tibble.
Exclude rows containing nan values.}
\usage{
turn_diff_in_tibble(Diff_loci)
}
\arguments{
\item{Diff_loci}{A list of matrices containing M, FST, G'ST and D values
for each loci in a sample. The values expected to be in the second column.}
}
\value{
A tibble object without nan values.
}
\description{
Transform a list of matrices containing differentiation statistics into a tibble.
Exclude rows containing nan values.
}
\examples{
library(tibble)
# Create example list of 2 matrices with 4 rows (M, FST, G'ST, D)
Diff_loci <- list(
  matrix(c("M", 0.85, NA, NA,
           "FST", 0.01, NA, NA,
           "G'ST", 0.02, NA, NA,
           "D", 0.015, NA, NA),
         ncol = 4, byrow = TRUE,
         dimnames = list(NULL, c("statistic", "value", "bound", "bound_closeness"))),

  matrix(c("M", 0.5, NA, NA,
           "FST", 0.02, NA, NA,
           "G'ST", 0.03, NA, NA,
           "D", 0.025, NA, NA),
         ncol = 4, byrow = TRUE,
         dimnames = list(NULL, c("statistic", "value", "bound", "bound_closeness")))
)

# Convert to tibble and remove rows with NaNs
clean_tibble <- turn_diff_in_tibble(Diff_loci)
print(clean_tibble)
}

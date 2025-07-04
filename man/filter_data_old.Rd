% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_subfraction_data.R
\name{filter_data_old}
\alias{filter_data_old}
\title{Filters out loci from other samples, duplicated rows, Clonal loci and loci with 0 in both subpopulations}
\usage{
filter_data_old(df, K = 2, filter_fixed = T)
}
\arguments{
\item{df}{A dataframe with loci in rows and different metrics in columns}

\item{K}{Number of subpopulations, K=2 by default}

\item{filter_fixed}{bool value, defines whether to filter out loci with 0 frequencies in both subpopulations}
}
\value{
A cleaned dataframe for further pipeline
}
\description{
Filters out loci from other samples, duplicated rows, Clonal loci and loci with 0 in both subpopulations
}
\examples{
library(dplyr)
# Simulated example dataset
df <- tibble::tibble(
  chr = c("1", "1", "1"),
  pos = c(12345, 12346, 12347),
  sub.frac.Sample1 = c(0.2, 0, 0.5),
  sub.frac.Sample2 = c(0.3, 0, 0.5),
  Clonal = c(FALSE, TRUE, FALSE),
  Sample = c("Sample1", "Sample1", "Sample2")
)
filtered_data <- filter_data_old(df, K=2)
}

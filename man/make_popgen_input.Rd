% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_popgen_input.R
\name{make_popgen_input}
\alias{make_popgen_input}
\title{Produces a list of allele frequency matrices for each locus,
formatted for input to the Diff() function}
\usage{
make_popgen_input(data_frequencies)
}
\arguments{
\item{data_frequencies}{A data frame of allele frequencies with dimensions (L × K),
where L is the number of loci (rows), and K is the number of subpopulations (columns)}
}
\value{
A list of K × 2 matrices of allele frequencies, one per locus.
Each matrix contains frequencies of the reference and alternative alleles
for each subpopulation.
}
\description{
Produces a list of allele frequency matrices for each locus,
formatted for input to the Diff() function
}
\examples{
# Create a toy dataframe of allele frequencies for 3 loci and 2 subpopulations
df <- data.frame(
  Subpop1 = c(0.8, 0.6, 0.4),
  Subpop2 = c(0.7, 0.5, 0.3)
)
# Use the function to generate a list of 2x2 frequency matrices
freq_list <- make_popgen_input(df)
}

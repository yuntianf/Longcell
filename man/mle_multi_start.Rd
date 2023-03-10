% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLRT_functions.R
\name{mle_multi_start}
\alias{mle_multi_start}
\title{mle for beta-binomial}
\usage{
mle_multi_start(exon_count, gene_count, start)
}
\arguments{
\item{exon_count}{a vector of exon count}

\item{gene_count}{a vector of gene count with the same length as exon_count}

\item{start}{a list of start points, each line should be a pair of alpha and beta}
}
\value{
a list of result for mle
}
\description{
mle for beta-binomial with multiple start points to avoid local optimization
}

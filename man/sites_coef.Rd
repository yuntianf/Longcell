% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splice_table.R
\name{sites_coef}
\alias{sites_coef}
\title{sites_coef}
\usage{
sites_coef(a, b, count)
}
\arguments{
\item{a}{A numerical vector for the existence of an exon in each isoform}

\item{b}{A numerical vector for the existence of an exon in each isoform}

\item{count}{The read count}
}
\value{
TRUE or FALSE, indicating if the two exons should be merged
}
\description{
judge if two sites can be merged
}
\details{
If two exons are co-existing or mutually exclusive, they have the same information,
and can considered as one exon in downstream analysis
}

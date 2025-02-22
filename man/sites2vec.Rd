% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splice_table.R
\name{sites2vec}
\alias{sites2vec}
\title{sites2vec}
\usage{
sites2vec(start, in_site, out_site, end, polyA, sites_index, strand)
}
\arguments{
\item{start}{A number recording the start position for a read}

\item{in_site}{A string vector recording the start position for each exon for a read}

\item{out_site}{A string vector recording the end position for each exon for a read}

\item{end}{A number recording the end position for a read}

\item{polyA}{A bool variable to indicate if the read has a polyA}

\item{sites_index}{A string vector recording all sites across all reads in this gene}

\item{strand}{The strand of this gene}
}
\value{
A binary vector indicating the splicing sites existence for a read
}
\description{
convert the splicing sites vector to a binary vector
}
\details{
convert the splicing sites vector to a binary vector, the spliced-in sites will be set as 1,
the spliced-out sites will be set as 0, the depleted sites which may be due to technical truncations
are set as NA.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exon_table.R
\name{extractExonTable}
\alias{extractExonTable}
\title{build an exon count table for a gene}
\usage{
extractExonTable(
  spliceOb,
  gene,
  gene_bed = NULL,
  gtf = NULL,
  cells = "all",
  exons = "all",
  exon_len_thresh = 10,
  sep = "|"
)
}
\arguments{
\item{spliceOb}{the input splice object}

\item{gene}{name of the target gene}

\item{gene_bed}{gene bed annotation}

\item{gtf}{gtf annotation}

\item{cells}{cells to use in building exon count table}

\item{exons}{exon to use in building exon count table}

\item{exon_len_thresh}{the minimum length of an exon to be preserved in exon count table}

\item{sep}{The seperation character between exons in the isoform string}
}
\value{
an exon count table for isoforms in each cell
}
\description{
build an exon count table for a gene, each row is an isoform in a cell
and each column is a meta-exon
}
\details{
build an exon count table for a gene, each row is an isoform in a cell
and each column is a meta-exon
}

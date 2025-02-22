% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phi_visual.R
\name{phiPlot}
\alias{phiPlot}
\title{phiPlot}
\usage{
phiPlot(
  data,
  psi_range = c(0.3, 0.7),
  phi_conf_upr = 0.2,
  psi_col = "mean_psi",
  phi_col = "phi",
  phi_conf_col = "phi_conf",
  annot_col = "symbol",
  annot_size = 4
)
}
\arguments{
\item{data}{the dataframe of phi, which should include columns for mean psi, phi, phi confidence
interval and the exon name}

\item{psi_range}{Only exons with mean psi within this range will be annotated}

\item{phi_conf_upr}{Only exons with phi confidence interval lower than this threshold
will be annotated}

\item{psi_col, phi_col, phi_conf_col, annot_col}{The colnames to indicate the mean psi,phi,phi confidence
interval and exon annotation within the input data}

\item{annot_size}{the text size for exon annotation}
}
\description{
plot mean psi vs. phi scatter plot for each exon, colored by the confidence
interval of phi
}

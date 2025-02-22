% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLRT_functions.R
\name{GLRT_sig_plot}
\alias{GLRT_sig_plot}
\title{plot for GLRT significant results}
\usage{
GLRT_sig_plot(
  data,
  gene_col = "gene",
  color_col = NULL,
  mean_diff = "mean_diff",
  var_diff = "var_diff",
  q = "q",
  mean_diff_thresh = 0.15,
  var_diff_thresh = 0.1,
  log_q_thresh = 3
)
}
\arguments{
\item{data}{the GLRT result, a dataframe output from genes_groups_GLRT()}

\item{gene_col}{The name of column in the data which stores the gene name}

\item{mean_diff}{The name of column in the data which stores the mean difference
of the psi distributions}

\item{var_diff}{The name of column in the data which stores the variance difference
of the psi distributions}

\item{q}{The name of column in the data which stores the q value}

\item{mean_diff_thresh}{the minimum of mean difference of psi to show the gene name in the scatter plot}

\item{var_diff_thresh}{the minimum of variance difference of psi to show the gene name in the scatter plot}

\item{q_thresh}{the minimum of -log(q,10) to show the gene name in the scatter plot}
}
\value{
a ggplot object
}
\description{
plot the scatter plot for mean change versus variance change of exon psi
between two cell groups
}

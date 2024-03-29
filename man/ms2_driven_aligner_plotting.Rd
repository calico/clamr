% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignment_plots.R
\name{ms2_driven_aligner_plotting}
\alias{ms2_driven_aligner_plotting}
\title{MS2 Driven Aligner Plotting}
\usage{
ms2_driven_aligner_plotting(
  rt_plot_data,
  print_plots,
  show_specific_samples = NULL
)
}
\arguments{
\item{rt_plot_data}{plot_data produced using \code{\link{ms2_driven_aligner}}.}

\item{print_plots}{If FALSE then return individual plots in a list; if TRUE then print a plot and silently return 1.}

\item{show_specific_samples}{A vector of sampleIds or if NULL a set of sampleIds will be randomely sampled. If sampleIds are named then they will be used for plots.}
}
\description{
MS2 Driven Aligner Plotting
}
\details{
\itemize{
    \item{1: dendrogram comparing samples based on shared MS2 manhattan distance}
    \item{2: tile plot comparing samples based on shared MS2 counts}
    \item{3: visualize RT deviation trend for individual samples}
    \item{4: compare the RT of MS/MS2 matched samples for pairs of samples}
    \item{5: consensus ppm drift of shown samples}
    \item{6: mass accuracy deviations of individual samples}
    }
}

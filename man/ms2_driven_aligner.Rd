% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignment.R
\name{ms2_driven_aligner}
\alias{ms2_driven_aligner}
\title{MS2 Driven Aligner}
\usage{
ms2_driven_aligner(
  mzroll_db_con,
  clamr_config,
  group_by_charge = FALSE,
  peak_quality_cutoff = 0.5,
  maximum_mz_set_size = 1000,
  cosine_cutoff = 0.95,
  sd_rt_resid_cutoff = 0.1,
  spline_ridge_penalty = 200,
  spline_degree = 4L,
  return_plot_data = FALSE,
  quietly = FALSE
)
}
\arguments{
\item{mzroll_db_con}{a connection to a mzroll database as produced by \code{\link{mzroll_db_sqlite}}}

\item{clamr_config}{a named list of mass spec parameters with special formatting of instrument tolerances generated by \code{\link{build_clamr_config}}.}

\item{group_by_charge}{Require that peaks match in precursor charge in order to be a possible match -- discard all peaks with unknown charge.}

\item{peak_quality_cutoff}{Minimum quality for a scan's matching peak for the ion to be used for alignment.}

\item{maximum_mz_set_size}{The maximum number of possibly matching ms2 events that will be considered.}

\item{cosine_cutoff}{Minimum cosine similarity between a pair of MS2 fragment profiles to group them into a common cluster.}

\item{sd_rt_resid_cutoff}{Cutoff for excluding a compound based on the standard deviation of its residuals (rt - fitted rt): sd(resid)/range(rt) < \code{sd_rt_resid_cutoff}.}

\item{spline_ridge_penalty}{Ridge penalty on sample-level splines which model deviations between observed retention times and the consensus rt of a compound. Used for the H parameter in \link[mgcv]{gam}.}

\item{spline_degree}{Degree of spline used to estimate drift (integer).}

\item{return_plot_data}{return data which can be used to plot alignment summaries with \code{\link{ms2_driven_aligner_plotting}}}

\item{quietly}{Hide messages and warnings}
}
\value{
an update mzroll_db_con where all features have had a function applied of the form RT_updated = RT_original + g(RT_original), where the aim of the g(RT_original) is to estimate retention time deviations between individual samples and a consensus sample.
}
\description{
Across a dataset, use MS1 and MS2 data to generate update function which can be used to transform retention times and mass accuracy in order to improve consistency.
}
\details{
This function uses MS1 and MS2 similarity to group MS2 events from different samples together and then determines the extent to which
individual samples systematically deviate in either retention time or mass accuracy.
}
\examples{
\dontrun{
# TO DO - add a smaller dataset so this would run fast enough for testing
library(dplyr)

mzroll_db_con <- clamshell::get_clamr_assets("mzkit_peakdetector.mzrollDB") \%>\%
  mzroll_db_sqlite()

clamr_config <- build_clamr_config(list(MS1tol = "10ppm", MS2tol = "20ppm"))
ms2_driven_aligner(mzroll_db_con, clamr_config)
}

}

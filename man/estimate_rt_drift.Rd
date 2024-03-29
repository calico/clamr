% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignment.R
\name{estimate_rt_drift}
\alias{estimate_rt_drift}
\title{Estimate RT Drift}
\usage{
estimate_rt_drift(
  sample_matched_compounds,
  previous_rt_update_key = NULL,
  sd_rt_resid_cutoff = 0.1,
  spline_ridge_penalty = 200,
  spline_degree = 4L,
  return_plot_data = FALSE
)
}
\arguments{
\item{sample_matched_compounds}{tibble containing at least sampleId, compound_id and rt.}

\item{previous_rt_update_key}{previous mapping of original RTs in raw files to RTs used in peaks and scans.}

\item{sd_rt_resid_cutoff}{Cutoff for excluding a compound based on the standard deviation of its residuals (rt - fitted rt): sd(resid)/range(rt) < \code{sd_rt_resid_cutoff}.}

\item{spline_ridge_penalty}{Ridge penalty on sample-level splines which model deviations between observed retention times and the consensus rt of a compound. Used for the H parameter in \link[mgcv]{gam}.}

\item{spline_degree}{Degree of spline used to estimate drift (integer).}

\item{return_plot_data}{return data which can be used to plot alignment summaries with \code{\link{ms2_driven_aligner_plotting}}}
}
\value{
a list containing sample specific gam fits and a sample independent fit which should be used to correct sample-level overfitting
}
\description{
Generate a mapping of observed sample retention times onto imputed (consistent) retention time.
}
\details{
Alternate between estimating compound-specific retention time and sample-specific retention time deviations.
}

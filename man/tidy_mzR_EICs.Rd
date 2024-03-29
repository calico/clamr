% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy_mzR_utils.R
\name{tidy_mzR_EICs}
\alias{tidy_mzR_EICs}
\title{Generated Extracted Ion Chromatographs using mzR}
\usage{
tidy_mzR_EICs(ms_file_path, mzroll_db_path = NULL, mz_tbl, clamr_config)
}
\arguments{
\item{ms_file_path}{Path to an mzXML, mzML, netCDF or mzData file}

\item{mzroll_db_path}{path to a .mzrollDB SQLite file.}

\item{mz_tbl}{a tibble/data.frame containing one row per m/z of interest, with this mz being specified by a variable "mz"}

\item{clamr_config}{a named list of mass spec parameters with special formatting of instrument tolerances generated by \code{\link{build_clamr_config}}.}
}
\value{
mz_tbl with a nested list of EICs named "eic"
}
\description{
Extract ion abundances per scan for a set of m/zs from a single sample.
}

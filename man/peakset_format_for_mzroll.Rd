% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peaksets.R
\name{peakset_format_for_mzroll}
\alias{peakset_format_for_mzroll}
\title{Peakset Format for Mzroll}
\usage{
peakset_format_for_mzroll(peakset)
}
\arguments{
\item{peakset}{generated by \code{\link{extract_peakset}}}
}
\value{
a peakset
}
\description{
Take a peakset object and remove fields which are not part of the core mzrollDB schema
and ensure that all required fields are present
}

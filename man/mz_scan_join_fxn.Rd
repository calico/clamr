% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/join_searches.R
\name{mz_scan_join_fxn}
\alias{mz_scan_join_fxn}
\title{Match scans to peaks for a single sample}
\usage{
mz_scan_join_fxn(sample_peaks, sample_scans)
}
\arguments{
\item{sample_peaks}{peaks for a single sample
columns: peakId, groupId, rt, rtmin, rtmax, peakMz, mzmin, mzmax, minscan, maxscan, peakAreaTop, quality}

\item{sample_scans}{scans for a single sample
columns: scan, rt, precursorMz, precursorCharge, precursorIC, precursorPurity, data}
}
\value{
scans matched to the best peak
}
\description{
Matching peaks [mzmin - mzmax, scanmin - scanmax] to a scan [mz, scan]
}

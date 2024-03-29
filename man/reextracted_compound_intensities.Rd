% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{reextracted_compound_intensities}
\alias{reextracted_compound_intensities}
\title{Reextracted Peak Intensity Set}
\usage{
reextracted_compound_intensities(
  curated_mzrolldb,
  full_mzrolldb,
  sample_folder,
  min_mz_delta = 0.01,
  min_rt_delta = 0.1,
  verbose = TRUE,
  debug = FALSE
)
}
\arguments{
\item{curated_mzrolldb}{mzrollDB file containing manually verified peakgroups.
peaks and peakgroups from this table are used in this method.}

\item{full_mzrolldb}{mzrollDB file containing original rt_update_key RT alignment table.
This table is used for transforming sample_rt <--> aligned_rt}

\item{sample_folder}{folder containing all raw samples to apply reextraction.}

\item{min_mz_delta}{adjust min and max m/z so that the difference is at least this value.}

\item{min_rt_delta}{adjust min and max rt so that the difference is at least this value.}

\item{verbose}{flag for explicit printing}

\item{debug}{flag for debugging}
}
\value{
intensities_table

   Each row indicates quant values for a single peak (an annotated feature in a single sample.)
   Both aligned RT coordinates and sample RT coordinates are given.

   \item{sample: }{sample name}
   \item{sampleId: }{ID number for sample in mzrolldb files}
   \item{groupId: }{group ID number from curated_mzrolldb file}
   \item{compoundName: }{String representation of annotated feature}
   \item{intensity: }{curated_mzrollDB-recorded intensity (peakAreaTop), or reextracted intensity}
   \item{is_identified: }{if TRUE, using quant from curated_mzrolldb. Else, reextract from sample file}
   \item{mzmin: }{minimum m/z associated with peak}
   \item{mzmax: }{maximum m/z associated with peak}
   \item{aligned_rtmin: }{minimum RT value associated with peak, in aligned RT space. Relies on rt_update_key table in full_mzrolldb}
   \item{aligned_rtmax: }{maximum RT value associated with peak, in aligned RT space. Relies on rt_update_key table in full_mzrolldb}
   \item{sample_rtmin: }{minimum RT value associated with peak, in mzML/mzXML sample space.}
   \item{sample_rtmax: }{maximum RT value associated with peak, in mzML/mzXML sample space.}
}
\description{
Produce a table of reextracted peak intensities based on a short
list of curated peak groups
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{strip_ms_file_extension_from_extensions}
\alias{strip_ms_file_extension_from_extensions}
\title{package-local function that does the actual work of modifying file extensions.}
\usage{
strip_ms_file_extension_from_extensions(file_names, file_extensions)
}
\arguments{
\item{file_names}{list of file names. can be a mixture of types (mzML, mzXML, etc)}

\item{file_extensions}{a character vector of extensions as generated by \link{get_valid_ms_file_extensions}.}
}
\value{
file with file extension stripped off.
}
\description{
package-local function that does the actual work of modifying file extensions.
}

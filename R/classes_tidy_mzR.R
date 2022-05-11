#' Read a mass spec file with mzR
#'
#' This function wraps and reformats \link[mzR]{openMSfile} converting its output to a list of tibbles.
#' It can open any of the file formats of openMSfile (mzXML, mzML, netCDF, mzData) as well as g-zipped versions of these formats.
#'
#' @details
#' If a \code{mzroll_db_path} is provided, then use this file's rt_update_key table (if it exists) to update retention times.
#'
#' @param ms_file_path Path to an mzXML, mzML, netCDF or mzData file
#' @inheritParams mzroll_db_sqlite
#' @param msLevels Filter for specific MS levels (if provided, this must be a vector of integers)
#'
#' @return tidy_mzR:
#' a list containing two tibbles:
#' \itemize{
#'   \item{header: summary of scan; one scan per row}
#'   \itemize{
#'     \item{acquisitionNum: scan #}
#'     \item{polarity: positive (1) or negative (0) mode}
#'     \item{other scan attributes}
#'   }
#'   \item{scan_data: abundances of ions in each scan; contains three variables: mz, ic, and scan}
#'   }
#'
#' @seealso \code{\link{tidy_mzR_plot_scans}}, \code{\link{tidy_mzR_plot_major_features}} and \code{\link{tidy_mzR_extract_mzs}}
#'
#' @export
tidy_mzR_from_msfile <- function(ms_file_path, mzroll_db_path = NULL, msLevels = NULL) {

  # test inputs
  stopifnot(class(ms_file_path) == "character")
  if (length(ms_file_path) != 1) {
    stop("ms_file_path must be a length 1 character vector")
  }

  if (!file.exists(ms_file_path)) {
    stop("No file with ms_file_path: ", ms_file_path)
  }

  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop('The "mzR" package must be installed to use this function',
      call. = FALSE
    )
  }

  valid_ms_file_extensions <- get_valid_ms_file_extensions()

  # test for valid extensions
  file_extension <- stringr::str_extract(ms_file_path, "\\.([a-zA-Z\\.]+)$")
  if (!(tolower(file_extension) %in% tolower(valid_ms_file_extensions))) {
    stop("The ms_file_path: ", ms_file_path, " does not have a valid extension\nvalid extensions are: ", paste(valid_ms_file_extensions, collapse = ", "))
  }

  stopifnot(
    length(class(msLevels)) >= 1,
    class(msLevels) %in% c("NULL", "integer")
  )

  a_file <- mzR::openMSfile(ms_file_path)

  scan_header <- mzR::header(a_file) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(scan = 1:dplyr::n())

  if (class(msLevels) != "NULL") {
    scan_header <- scan_header %>%
      dplyr::filter(msLevel %in% msLevels)

    if (nrow(scan_header) == 0) {
      warning("No scans with an msLevel of ", paste(msLevels, collapse = " or "))
      return(NULL)
    }
  }

  if (!any(c("NULL", "character") %in% class(mzroll_db_path))) {
    stop(mzroll_db_path, 'is of class "', class(mzroll_db_path), '" mzroll_db_path must either be a length 1 character path or NULL')
  }

  if (class(mzroll_db_path) == "character") {
    stopifnot(length(mzroll_db_path) == 1)

    if (!file.exists(mzroll_db_path)) {
      stop("mzroll_db_path not found at: ", mzroll_db_path)
    }

    mzroll_db_con <- mzroll_db_sqlite(mzroll_db_path)
    if ("rt_update_key" %in% DBI::dbListTables(mzroll_db_con)) {

      # find the sampleId corresponding to the ms_file_path sample
      ms_file_sampleId <- match_mzroll_db_to_file(mzroll_db_con, ms_file_path)

      if (class(ms_file_sampleId) == "NULL") {
        warning(ms_file_path, " does not match any sample names in ", mzroll_db_path, "; no retention time update can be performed")
      } else {

        # apply retention time update using rt_update_key for matched sampleId
        sample_rt_update_key <- dplyr::tbl(mzroll_db_con, "rt_update_key") %>%
          dplyr::filter(sampleId == ms_file_sampleId) %>%
          dplyr::collect()

        if (nrow(sample_rt_update_key) == 0) {
          warning("sample rt_udpate_key contains no information for sampleId: ", ms_file_sampleId)
        }

        # update scan header retention time using stored rt alignment
        scan_header$retentionTime <- mzkitcpp::update_rts(sample_rt_update_key, scan_header$retentionTime, rep(ms_file_sampleId, nrow(scan_header)), debug = F)$updated_rts
      }
    } else {
      warning('"rt_update_key" not present in mzrollDB; no retention time update will be performed')
    }
  }
  # do nothing if mzroll_db_path is NULL

  scan_data <- mzR::peaks(a_file, scans = scan_header$scan) %>%
    purrr::map2(scan_header$scan, cbind) %>%
    purrr::map(function(x) {
      colnames(x) <- c("mz", "ic", "scan")
      x
    }) %>%
    purrr::map(tibble::as_tibble) %>%
    dplyr::bind_rows()

  tidy_mzR <- list(
    header = scan_header,
    scan_data = scan_data
  )

  class(tidy_mzR) <- c("tidy_mzR", "list")

  tidy_mzR
}

tidy_mzR_test <- function(tidy_mzR) {
  if (!("tidy_mzR" %in% class(tidy_mzR))) {
    stop("\"tidy_mzR\" is not a tidy_mzR class; generate this object using tidy_mzR_from_msfile()")
  }

  if (!("list" %in% class(tidy_mzR))) {
    stop("\"tidy_mzR\" is malformed; generate this object using tidy_mzR_from_msfile()")
  }

  if (any(names(tidy_mzR) != c("header", "scan_data"))) {
    stop("\"tidy_mzR\" is malformed; generate this object using tidy_mzR_from_msfile()")
  }

  invisible(0)
}

summary.tidy_mzR <- function(tidy_mzR) {
  tidy_mzR_test(tidy_mzR)

  summary_attributes <- tibble::tribble(
    ~"Attribute", ~"Value",
    "Max RT", round(max(tidy_mzR$header$retentionTime / 60), 3),
    "Min MZ", round(max(tidy_mzR$header$lowMZ), 5),
    "Max MZ", round(max(tidy_mzR$header$highMZ), 5)
  )

  scan_attributes <- tidy_mzR$header %>%
    dplyr::count(msLevel) %>%
    dplyr::rename(Attribute = msLevel, Value = n) %>%
    dplyr::mutate(
      Attribute = as.character(glue::glue("MS {Attribute} scans")),
      Value = Value
    )

  do.call(dplyr::bind_rows, list(
    summary_attributes,
    scan_attributes
  )) %>%
    knitr::kable() %>%
    print()
}

#' Match .mzrollDB to file
#'
#' find the sampleId corresponding to a provided sample file path
#'
#' @inheritParams test_mzroll_db_con_schema
#' @inheritParams tidy_mzR_from_msfile
#'
#' @return a matched sampleId or NULL if no match was found
match_mzroll_db_to_file <- function(mzroll_db_con, ms_file_path) {
  if (!("samples" %in% DBI::dbListTables(mzroll_db_con))) {
    stop('"samples" table not present in mzrollDB; mzrollDB is mal-formed, samples should always be present')
  }

  samples <- dplyr::tbl(mzroll_db_con, "samples") %>%
    dplyr::select(sampleId, name) %>%
    dplyr::collect()

  stopifnot(class(ms_file_path) == "character", length(ms_file_path) == 1, file.exists(ms_file_path))

  ms_filename <- basename(ms_file_path)

  if (ms_filename %in% samples$name) {
    output <- samples %>%
      dplyr::filter(name == ms_filename) %>%
      {
        .$sampleId[1]
      }
  } else {
    output <- NULL
  }

  return(output)
}

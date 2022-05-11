#' Match scans to peaks for a single sample
#'
#' Matching peaks [mzmin - mzmax, scanmin - scanmax] to a scan [mz, scan]
#'
#' @param sample_peaks peaks for a single sample
#'      columns: peakId, groupId, rt, rtmin, rtmax, peakMz, mzmin, mzmax, minscan, maxscan, peakAreaTop, quality
#'
#' @param sample_scans scans for a single sample
#'      columns: scan, rt, precursorMz, precursorCharge, precursorIC, precursorPurity, data
#'
#' @return scans matched to the best peak
mz_scan_join_fxn <- function(sample_peaks, sample_scans) {

  # Handle null cases
  if (is.null(sample_scans)) {
    return(tibble::tibble(
      "scan" = integer(0),
      "rt" = numeric(0),
      "precursorMz" = numeric(0),
      "precursorCharge" = integer(0),
      "precursorIc" = numeric(0),
      "precursorPurity" = numeric(0),
      "data" = character(0)
    ))
  } else if (is.null(sample_peaks)) {
    return(sample_scans)
  }

  # join peaks and scans if minscan(peak) <= scan(scan) <= maxscan(peak)
  scan_match <- sample_scans %>%
    dplyr::select(scan, precursorMz) %>%
    # expand join window by 1 since the scan where an MS2 is produced will not also contain MS1 data
    dplyr::mutate(minscan = scan - 1, maxscan = scan + 1) %>%
    fuzzyjoin::interval_inner_join(
      sample_peaks %>%
        dplyr::select(peakId, peakMz, mzmin, mzmax, minscan, maxscan),
      by = c("minscan", "maxscan"), maxgap = 0
    )

  # filter for mz match if mzmin(peak) <= precursorMz(scan) <= mzmax(peak)
  best_full_match <- scan_match %>%
    dplyr::filter(precursorMz >= mzmin - 0.001, precursorMz <= mzmax + 0.001) %>%
    dplyr::group_by(scan) %>%
    # dplyr::filter(dplyr::n() > 1) %>%
    dplyr::arrange(abs(precursorMz - peakMz)) %>%
    dplyr::slice(1) %>%
    dplyr::select(scan, peakId)

  sample_scans %>%
    dplyr::left_join(best_full_match, by = "scan") %>%
    dplyr::left_join(sample_peaks, by = "peakId")
}

#' Match peaks to peaks for a single sample
#'
#' @param sample_peaks1 peaks for a single sample (mzmin, mzmax, minscan, maxscan)
#' @param sample_peaks2 peaks for a single sample (mzmin, mzmax, minscan, maxscan)
#'
#' Matching peaks [mzmin - mzmax, minscan - maxscan] to peaks [mzmin - mzmax, minscan - maxscan]
mz_peak_join_fxn <- function(sample_peaks1, sample_peaks2) {
  # join peaks if scan overlap exists
  sample_peaks1 %>%
    fuzzyjoin::interval_inner_join(sample_peaks2, by = c("minscan", "maxscan"), maxgap = 0) %>%
    # filter for mz interval overlap
    dplyr::filter(mzmin.x <= mzmax.y & mzmin.y <= mzmax.x)
}

#' Join peaks to scans when their is a scan range and mz match for a given sample
#'
#' @inheritParams test_mzroll_db_con_schema
#'
#' @export
join_peaks_to_scans <- function(mzroll_db_con) {
  peak_data <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT peakId, groupId, sampleId, rt, rtmin, rtmax, peakMz, mzmin, mzmax, minscan, maxscan, peakAreaTop, quality FROM peaks")) %>%
    dplyr::collect() %>%
    dplyr::rename(peakRt = rt)

  scan_data <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT sampleId, scan, rt, precursorMz, precursorCharge, precursorIC, precursorPurity, data FROM scans WHERE mslevel == 2")) %>%
    dplyr::collect()

  peaks_joined_to_scans <- tidyr::nest(peak_data, peaks = -sampleId) %>%
    dplyr::left_join(tidyr::nest(scan_data, scans = -sampleId), by = "sampleId") %>%
    dplyr::mutate(matches = purrr::map2(peaks, scans, mz_scan_join_fxn)) %>%
    dplyr::select(-peaks, -scans) %>%
    tidyr::unnest(matches)

  peaks_joined_to_scans
}

#' Augment Peak Groups
#'
#' Extract peakgroups table and add peakgroup level summaries from peaks
#'
#' @inheritParams test_mzroll_db_con_schema
#'
#' @export
augment_peakgroups <- function(mzroll_db_con) {
  peakgroups <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT * FROM peakgroups")) %>%
    dplyr::collect()

  peak_summaries <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT groupId, rt as peakRt, peakMz, peakAreaTop, quality FROM peaks")) %>%
    dplyr::collect() %>%
    dplyr::group_by(groupId) %>%
    dplyr::summarize(
      mean_peakMz = sum(peakMz * peakAreaTop) / sum(peakAreaTop),
      mean_peakRt = sum(peakRt * peakAreaTop) / sum(peakAreaTop),
      mean_peakAreaTop = mean(peakAreaTop),
      median_peakAreaTop = stats::median(peakAreaTop),
      mean_peakQuality = mean(quality)
    )

  peakgroups %>%
    dplyr::left_join(peak_summaries, by = "groupId")
}

#' MS2 Scans from files
#'
#' Get all MS2 scans from all samples in a directory
#'
#' @param mzroll_db_con connection to mzroll database
#' @param sample_dir directory containing mzML / mzXML raw files
#'
#' @export
ms2_scans_from_raw_files <- function(mzroll_db_con, sample_dir = NULL) {

  # Use current working directory if no sample dir provided, or sample_dir argument not a string
  if (is.null(sample_dir) || class(sample_dir) != "character") {
    sample_dir <- "."
  }

  sample_data <- dplyr::tbl(mzroll_db_con, dplyr::sql("SELECT sampleId, name, filename FROM samples")) %>%
    dplyr::collect()

  samples <- list.files(sample_dir, pattern = "*.mzML|.mzXML", full.names = TRUE)

  full_path <- function(name) {
    coord <- grep(name, samples)
    if (length(coord) > 1) {
      stop("Cannot map sample names from directory to names stored in mzrollDB file!")
    }
    samples[coord]
  }

  sample_data <- sample_data %>%
    dplyr::mutate(filename = sapply(name, full_path, USE.NAMES = F))

  scan_data <- mzkitcpp::get_ms2_scans(sample_data, debug = F)

  # fields used from SQLite select statement
  scan_data <- scan_data %>%
    dplyr::select(sampleId, scan, rt, precursorMz, precursorCharge, precursorIc, precursorPurity, data)

  scan_data
}

#' Join GroupIds to {mz/z, rt}
#'
#' Find the most likely groupId corresponding to a set of mz, rt tuples.
#'
#' @param mz_rts a tibble containing mz and rt
#' @inheritParams test_mzroll_db_con_schema
#' @inheritParams test_clamr_config
#'
#' @return mz_rts with groupIds added
#'
#' @export
join_groupIds_to_mz_rts <- function(mz_rts, mzroll_db_con, clamr_config) {
  stopifnot("data.frame" %in% class(mz_rts), all(c("mz", "rt") %in% colnames(mz_rts)))
  clamr::require_tolerances(clamr_config, 1L)

  mzroll_peakgroups <- clamr::extract_peakset(mzroll_db_con)$peakgroups %>%
    dplyr::select(groupId, mz = mean_peakMz, rt = mean_peakRt)


  variable_tolerances <- tibble::tribble(
    ~variable, ~tolerance, ~relative_or_absolute,
    "mz", clamr_config$MS1tol$tol, clamr_config$MS1tol$absolute_or_relative,
    "rt", 1, "absolute"
  )

  matched_peakgroups <- clamr::join_matching_standards(mz_rts,
    mzroll_peakgroups,
    variable_tolerances = variable_tolerances,
    threshold = 1
  ) %>%
    dplyr::arrange(match_distance) %>%
    dplyr::group_by(query_mz, query_rt) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(mz = query_mz, rt = query_rt) %>%
    dplyr::select(!!!rlang::quos(c(colnames(mz_rts), "groupId")))

  return(matched_peakgroups)
}

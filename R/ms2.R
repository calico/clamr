#' Extract fragments only
#'
#' Read out the ms2 data exactly as it was captured from the instrument.
#' @param ms2scan a single ms2 scan, as produced by nested_peakgroup_features$scans[[x]]
#'
#' @export
extract_fragments <- function(ms2scan) {
  val <- strsplit(ms2scan$data, split = "\\[|,|\\]")[[1]]
  val <- val[val != ""]
  fragments <- tibble::tibble(
    ms2_mz = as.numeric(val[seq_along(val) %% 2 == 1]),
    ms2_ic = as.numeric(val[seq_along(val) %% 2 == 0])
  )
}

#' Sum Spectra
#'
#' @param spectra_set a tibble containing mz, ic, and optionally label.
#' @param MS1tol a tolerance list created with \code{format_mass_accuracy_input}.
#'
#' @return a tibble containing mz, ic, & optionally label
#'
#' @export
sum_spectra <- function(spectra_set, MS1tol) {
  checkmate::assertDataFrame(spectra_set)
  stopifnot(c("mz", "ic") %in% colnames(spectra_set))
  stopifnot(class(MS1tol) == "list", all(c("tol", "absolute_or_relative") %in% names(MS1tol)))

  if (!"label" %in% colnames(spectra_set)) {
    spectra_set <- spectra_set %>% dplyr::mutate(label = NA_character_)
  }

  summed_spectra <- spectra_set %>%
    dplyr::arrange(mz) %>%
    # bin mzs based on mass tolerance
    dplyr::mutate(mz_bin = clamr::find_mz_jumps(
      .$mz,
      MS1tol$tol,
      MS1tol$absolute_or_relative,
      collapse = FALSE
    )) %>%
    dplyr::group_by(mz_bin) %>%
    dplyr::summarize(
      mz = sum(mz * ic) / sum(ic),
      ic = sum(ic),
      label = paste(
        unique(label[!is.na(label) & label != ""]),
        collapse = ", "
      ),
      .groups = "drop"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-mz_bin)

  if (!"label" %in% colnames(spectra_set)) {
    summed_spectra <- summed_spectra %>% dplyr::select(-label)
  }

  return(summed_spectra)
}

#' Aggregate Scans
#'
#' Combine multiple scans into a consensus spectrum
#'
#' @param group_scans a tibble which minimally contains precursorMz, precursorPurity, quality (derived from the parent peak), and data.
#' @inheritParams extract_and_group_fragments
#' @param n_top_spectra_summed integer counts of maximum number of spectra to aggregate
#' @param quality_weights length 2 named vector with names "purity" and "quality" indicating the relative amount to weight by precursor purity (i.e., the amount of isolated signal matching the precursorMz) versus peak quality (i.e., good peak shapes).
aggregate_scans <- function(group_scans, MS2tol, n_top_spectra_summed = 3L, quality_weights = c("purity" = 2, "quality" = 1)) {
  stopifnot("integer" %in% class(n_top_spectra_summed), length(n_top_spectra_summed) == 1)
  stopifnot(length(quality_weights) == 2, all(names(quality_weights) == c("purity", "quality")), all(quality_weights > 0))
  stopifnot(all(c("tol", "absolute_or_relative") %in% names(MS2tol)))

  best_peakgroup_fragmentations <- group_scans %>%
    # for a given peakgroup, prioritize the best spectra to sum based on MS2 purity (precursorPurity) & peak quality (quality)
    dplyr::mutate(MS2_quality = (precursorPurity * quality_weights["purity"] + quality * quality_weights["quality"]) / sum(quality_weights)) %>%
    dplyr::arrange(desc(MS2_quality)) %>%
    # take up to the top n_top_spectra_summed best MS2s
    dplyr::slice(1:n_top_spectra_summed) %>%
    dplyr::select(precursorMz, precursorPurity, quality, MS2_quality, data)

  best_peakgroup_fragData <- best_peakgroup_fragmentations %>%
    # unpack fragmentation data (and bin fragments for a given peak group based on MS2 tolerance)
    dplyr::mutate(peak_num = 1:dplyr::n()) %>%
    tidyr::nest_legacy(.key = "frag") %>%
    dplyr::mutate(fragData = purrr::map(frag, extract_and_group_fragments, MS2tol = MS2tol)) %>%
    dplyr::select(-frag) %>%
    tidyr::unnest(fragData)

  precursorMz <- mean(best_peakgroup_fragmentations$precursorMz)

  fragmentationData_summary <- best_peakgroup_fragData %>%
    # sum the top n_top_spectra_summed spectra
    dplyr::group_by(mz_set) %>%
    dplyr::summarize(
      mz = sum(ms2_mz * ms2_ic) / sum(ms2_ic),
      ic = mean(c(ms2_ic, rep(0, nrow(best_peakgroup_fragmentations) - dplyr::n()))),
      mean_ms2_ic_frac = mean(c(ms2_ic_frac, rep(0, nrow(best_peakgroup_fragmentations) - dplyr::n()))),
      cv_ms2_ic_frac = stats::sd(c(ms2_ic_frac, rep(0, nrow(best_peakgroup_fragmentations) - dplyr::n()))) / mean_ms2_ic_frac,
      fractionSamples = dplyr::n() / nrow(best_peakgroup_fragmentations)
    ) %>%
    dplyr::mutate(is_precursor = ifelse(abs(mz - precursorMz) < 0.01, TRUE, FALSE))

  list(
    best_peakgroup_fragmentations = best_peakgroup_fragmentations,
    fragmentationData_summary = fragmentationData_summary
  )
}

#' Extract and group fragments
#'
#' For a set of MS2 spectra pull out all MS2 fragments and group them into sets by fragment m/z
#'
#' @inheritParams find_consistent_ms2_fingerprints
#' @param MS2tol a list of MS2 tol parameters produced with \code{\link{build_clamr_config}}
#' @param inconsistency_proportion fraction of signal which can exist between peaks to still allow them to be resolved
#'
#' @return a tibble of an aggregated MS2 spectra
extract_and_group_fragments <- function(mz_set, MS2tol, inconsistency_proportion = 0.005) {
  # extract and order fragments
  lapply(mz_set$peak_num, function(a_peak) {
    # pull out the MS2 m/z and IC values for each peak (by row) in a set
    # this is definitely inelegant
    val <- strsplit(mz_set$data[a_peak], split = "\\[|,|\\]")[[1]]
    val <- val[val != ""]
    tibble::tibble(peak_num = a_peak, ms2_mz = as.numeric(val[seq_along(val) %% 2 == 1]), ms2_ic = as.numeric(val[seq_along(val) %% 2 == 0]))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(ms2_mz) %>%
    # group fragments and obtain a consensus
    # identify jumps between blocks of m/z values
    dplyr::mutate(mz_set = find_mz_jumps(
      sorted_mzs = .$ms2_mz,
      mass_accuracy_tolerance = MS2tol$tol,
      absolute_or_relative = MS2tol$absolute_or_relative,
      n_max_gap_inconsistencies = as.integer(floor(nrow(.) * inconsistency_proportion)), collapse = FALSE
    )) %>%
    # generate consistency within each mz_set (for a given sample) by summing signal and adopting a consensus m/z
    # some peaks within the same m/z groups have MS2 signal at multiple (very similar) m/z, so add their IC
    dplyr::group_by(peak_num, mz_set) %>%
    dplyr::summarize(
      ms2_mz = stats::median(ms2_mz),
      ms2_ic = sum(ms2_ic)
    ) %>%
    # for each mz_set take a consensus m/z
    dplyr::group_by(mz_set) %>%
    dplyr::mutate(ms2_mz = stats::median(ms2_mz)) %>%
    dplyr::group_by(peak_num) %>%
    dplyr::mutate(ms2_ic_frac = ms2_ic / sum(ms2_ic)) %>%
    dplyr::ungroup()
}

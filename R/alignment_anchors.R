#' Find Candidate Anchor Points
#'
#' Identify peaks which are exclusively associated with a given m/z and thus
#' can be unambiguously identified across samples for retention time alignment.
#'
#' @inheritParams test_mzroll_db_con_schema
#' @inheritParams test_clamr_config
#' @param n_anchors number of anchor points to search for.
#' @param anchor_rank_cutoff rank cutoff for a viable anchors, smaller is better, ranks are scaled to [0,1].
#'
#' @return a tibble containing the 10-best scoring {mz, rt} pairs within each rt bin.
#'
#' @details
#' The function aims to identify {mz, rt} pairs which can serve as anchor points
#' for retention time alignment. MAVEN will find these anchors in each sample
#' and enforce that all anchors peaks have the same retention time, other peaks
#' retention times will then be shifted by linearly interpolating between
#' flanking anchors.
#'
#' @examples
#' \dontrun{
#' future::plan("multicore")
#' mzroll_db_con <- clamr::mzroll_db_sqlite("/tmp/M005A-set123.mzrollDB")
#' clamr_config <- clamr::build_clamr_config(list(MS1tol = "10ppm", MS2tol = "20ppm"))
#' find_candidate_anchors(mzroll_db_con, clamr_config, n_anchors = 30L, anchor_rank_cutoff = 0.5)
#' }
#'
#' @export
find_candidate_anchors <- function(mzroll_db_con, clamr_config, n_anchors = 30L, anchor_rank_cutoff = 0.5) {
  require_tolerances(clamr_config, 1L)

  stopifnot(
    any(c("numeric", "integer") %in% class(n_anchors)),
    length(n_anchors) == 1,
    n_anchors > 0
  )
  if (class(n_anchors) == "numeric") {
    n_anchors <- as.integer(n_anchors)
  }

  stopifnot(
    "numeric" %in% class(anchor_rank_cutoff),
    length(anchor_rank_cutoff) == 1,
    anchor_rank_cutoff > 0,
    anchor_rank_cutoff < 1
  )

  # work with just peaks (so we aren't assuming good peak groups)

  peaks <- dplyr::tbl(mzroll_db_con, "peaks") %>%
    dplyr::collect()

  # find peaks within mass tolerance of query peak separately for each sample
  # generates two summaries - (1) the fraction of ic in the query peak relative
  # to all peaks in mass tolerance, and (2) the standard deviation of rt across
  # the peaks within mass tolerance as a measure of how strongly multiple
  # peaks are clustered together.

  peak_exclusivity <- peaks %>%
    tidyr::nest(sample_peaks = -sampleId) %>%
    dplyr::mutate(peak_ic_fractions = furrr::future_map(sample_peaks,
      calculate_peak_exclusivity,
      MS1tol = clamr_config$MS1tol
    )) %>%
    dplyr::select(peak_ic_fractions) %>%
    tidyr::unnest(peak_ic_fractions)

  peak_attributes <- peaks %>%
    # find all measures of peak quality and exclusivity
    dplyr::select(sampleId, peakId, peakMz, rt, peakAreaTop, quality) %>%
    dplyr::left_join(peak_exclusivity, by = "peakId") %>%
    # remove clearly invalid peaks due to abysmal quality or exclusivity
    dplyr::filter(
      weighted_rt_sd < 2,
      quality > 0.4,
      peakAreaTop > 2^12
    ) %>%
    # separate peaks into sets which can be trivially separated by m/z
    dplyr::arrange(peakMz) %>%
    dplyr::mutate(mz_set = find_mz_jumps(.$peakMz,
      clamr_config$MS1tol$tol,
      clamr_config$MS1tol$absolute_or_relative,
      collapse = FALSE
    )) %>%
    dplyr::ungroup()

  # remove mz_sets which include too many different masses (mean_total_fractions)
  # or too few distinct samples

  n_samples <- length(unique(peaks$sampleId))

  mz_set_coverage <- peak_attributes %>%
    dplyr::group_by(mz_set, sampleId) %>%
    dplyr::summarize(total_fractions = sum(peak_ic_fraction)) %>%
    dplyr::group_by(mz_set) %>%
    dplyr::summarize(
      max_total_fractions = max(total_fractions),
      n = dplyr::n()
    ) %>%
    dplyr::filter(
      max_total_fractions == 1,
      n > (n_samples * 0.8)
    )

  # for viable mz_sets, generate mz-level summaries of exclusivity

  exclusive_mzs <- peak_attributes %>%
    tidyr::nest(mz_set_data = -mz_set) %>%
    # drop mz_sets with few samples included
    dplyr::semi_join(mz_set_coverage, by = "mz_set") %>%
    dplyr::mutate(mz_exclusivity = purrr::map(mz_set_data, summarize_mz_set)) %>%
    dplyr::select(mz_set, mz_exclusivity) %>%
    tidyr::unnest(mz_exclusivity)

  # rank mz_sets for each measure so that low ranks are better and then
  # sum ranks to generate an overall exclusivity measure.

  # create a ranking from each measure
  exclusive_mz_ranks <- exclusive_mzs %>%
    dplyr::mutate(
      quality_h_mean_rank = percent_rank(desc(quality_h_mean)),
      peakAreaTop_h_mean_rank = percent_rank(desc(peakAreaTop_h_mean)),
      ic_fraction_h_mean_rank = percent_rank(desc(ic_fraction_h_mean)),
      rt_sd_mean_rank = percent_rank(rt_sd_mean),
      # up-weight overall representation
      n_rank = percent_rank(desc(n)) * 3,
      # summed and re-scaled overall rank
      overall_quality_rank = (quality_h_mean_rank + peakAreaTop_h_mean_rank +
        ic_fraction_h_mean_rank + rt_sd_mean_rank + n_rank) / 7
    )

  # find the best set of anchor points

  rt_breaks <- seq(0,
    max(peaks$rtmax),
    length.out = n_anchors + 1
  )

  candidate_anchors <- exclusive_mz_ranks %>%
    # filter based on user-specified rank cutoff
    dplyr::filter(overall_quality_rank <= anchor_rank_cutoff) %>%
    # find the rt bin that each peak falls within
    dplyr::mutate(rt_bin = cut(rt,
      breaks = rt_breaks,
      include.lowest = TRUE,
      right = FALSE
    ))

  message(glue::glue("{nrow(candidate_anchors)} candidate anchor point peaks were found
                     across {length(unique(candidate_anchors$rt_bin))} of the {n_anchors} target retention times"))

  top10_candidates <- candidate_anchors %>%
    dplyr::group_by(rt_bin) %>%
    dplyr::arrange(overall_quality_rank) %>%
    dplyr::slice(1:10) %>%
    dplyr::ungroup()

  return(top10_candidates)
}

calculate_peak_exclusivity <- function(sample_peaks, MS1tol) {
  variable_tolerances <- tibble::tribble(
    ~variable, ~tolerance, ~relative_or_absolute,
    "peakMz", MS1tol$tol, MS1tol$absolute_or_relative
  )

  # left join peaks on peaks based based on mass tolerance to find a set of
  # peakId_2s with similar mass to each peakId_1s.

  mz_matches <- join_matching_standards(
    sample_peaks %>%
      dplyr::select(
        peakId_1 = peakId,
        peakMz = peakMz,
        peakAreaTop_1 = peakAreaTop
      ),
    sample_peaks %>%
      dplyr::select(
        peakId_2 = peakId,
        peakMz = peakMz,
        peakAreaTop_2 = peakAreaTop,
        peakRt_2 = rt
      ),
    variable_tolerances
  )

  # for each peakId_1, find the total ic of all peakId_2s within mass tolerance.
  # also calculate the sd of peak rts weighted by ic

  peak_exclusivity <- mz_matches %>%
    dplyr::group_by(peakId_1) %>%
    dplyr::summarize(
      summed_ic_within_tol = sum(peakAreaTop_2),
      rt_w_mean = weighted.mean(peakRt_2, peakAreaTop_2),
      weighted_rt_sd = sqrt(sum(peakAreaTop_2 * (peakRt_2 - rt_w_mean)^2) / sum(peakAreaTop_2))
    )

  # find ic fraction of each peakId_1 relative to all peaks within tolerance
  # (including itself)

  peak_exclusivity <- sample_peaks %>%
    dplyr::left_join(peak_exclusivity, by = c("peakId" = "peakId_1")) %>%
    dplyr::mutate(peak_ic_fraction = peakAreaTop / summed_ic_within_tol) %>%
    dplyr::select(peakId, peak_ic_fraction, weighted_rt_sd)

  return(peak_exclusivity)
}

summarize_mz_set <- function(mz_set_data) {

  # take the highest signal peak within each mz_set
  mz_set_data %>%
    dplyr::arrange(dplyr::desc(peak_ic_fraction)) %>%
    dplyr::group_by(sampleId) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(
      peakMz = mean(peakMz),
      rt = mean(rt),
      quality_h_mean = 1 / mean(1 / quality),
      peakAreaTop_h_mean = 1 / mean(1 / peakAreaTop),
      ic_fraction_h_mean = 1 / mean(1 / peak_ic_fraction),
      rt_sd_mean = mean(weighted_rt_sd),
      n = dplyr::n()
    )
}

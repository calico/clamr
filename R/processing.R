#' Remove Redundant Peakgroups
#'
#' Find and remove redundant peakgroups based on MS1 concordance Mz, Rt and Quants
#'
#' @inheritParams test_mzroll_db_con_schema
#' @inheritParams test_clamr_config
#' @param rt_grouping_tol maximum departure between retention time of two groups to be consistent
#' @param group_corr_cutoff minimum disagreement in terms of correlation between two groups to be consistent
#'
#' @export
remove_redundant_peakgroups <- function(mzroll_db_con, clamr_config, rt_grouping_tol = 0.5, group_corr_cutoff = 0.98) {
  debugr::dwatch(msg = "started remove_redundant_peakgroups. [clamr<processing.R>::remove_redundant_peakgroups]\n")

  require_tolerances(clamr_config, 1L)
  stopifnot(class(rt_grouping_tol) == "numeric", length(rt_grouping_tol) == 1, rt_grouping_tol >= 0)
  stopifnot(class(group_corr_cutoff) == "numeric", length(group_corr_cutoff) == 1, group_corr_cutoff >= 0, group_corr_cutoff <= 1)

  # extract peakset: peaks and peakgroups
  peakset <- extract_peakset(mzroll_db_con)

  debugr::dwatch(msg = "created peakgroups_summarized. [clamr<processing.R>::remove_redundant_peakgroups]\n")

  # separate in groups which cannot be trivially separated by mass
  peakgroup_mzgroups <- peakset$peakgroups %>%
    dplyr::arrange(mean_peakMz) %>%
    dplyr::mutate(mz_set = find_mz_jumps(mean_peakMz,
      clamr_config$MS1tol$tol,
      clamr_config$MS1tol$absolute_or_relative,
      collapse = FALSE
    ))

  debugr::dwatch(msg = "created peakgroups_mzgroups. [clamr<processing.R>::remove_redundant_peakgroups]\n")

  compare_groups <- peakgroup_mzgroups %>%
    tidyr::nest(mz_set_data = -mz_set) %>%
    dplyr::mutate(mz_set_groups = purrr::map(mz_set_data, find_unique_peakgroups,
      peaks = peakset$peaks, rt_grouping_tol = rt_grouping_tol, group_corr_cutoff = group_corr_cutoff
    ))

  debugr::dwatch(msg = "created compare_groups. [clamr<processing.R>::remove_redundant_peakgroups]\n")

  unique_groupIds <- compare_groups %>%
    dplyr::select(-mz_set_data) %>%
    tidyr::unnest(mz_set_groups) %>%
    # take the highest quality group from each cluster
    dplyr::left_join(peakset$peakgroups %>% dplyr::select(groupId, mean_quality), by = "groupId") %>%
    dplyr::group_by(mz_set, cluster) %>%
    dplyr::arrange(desc(mean_quality)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  debugr::dwatch(msg = "created unique_group_Ids. [clamr<processing.R>::remove_redundant_peakgroups]\n")

  peakset$peakgroups %>%
    dplyr::anti_join(unique_groupIds, by = "groupId") %>%
    {
      cat(nrow(.), "redundant peakgroups removed\n")
    }

  # remove redundant groupId from mzroll_db_con
  sqlite_write(mzroll_db_con, "peakgroups", peakset$peakgroups %>% dplyr::semi_join(unique_groupIds, by = "groupId"), overwrite = TRUE)
  sqlite_write(mzroll_db_con, "peaks", peakset$peaks %>% dplyr::semi_join(unique_groupIds, by = "groupId"), overwrite = TRUE)

  debugr::dwatch(msg = "updated mzroll_db. [clamr<processing.R>::remove_redundant_peakgroups]\n")

  invisible(0)
}

#' Identify Split Peaks
#'
#' Find split peaks where a single analyte is distributed across one or
#' more peakgroups.
#'
#' @inheritParams test_mzroll_db_con_schema
#' @inheritParams test_clamr_config
#' @param max_rt_deviation maximum rt deviation between two peakgroups which
#'   could be the same analyte.
#' @param ic_floor floor all low or missing signals to this value for the
#'   purpose of calculating anti-correlation in log-space.
#' @param anticorrelation_co cutoff for what is a strong anti-correlation as a
#'   sign of peak splitting (most signal of some samples in group A and others in group B will induce a negative correlation between A and B).
#' @param signal_frac_IQR_co cutoff for interquartile range of fractional
#'   signal spread between candidate split peak pairs.
#' @param spectra_corr_co cutoff for correlation between consensus spectra
#'   of the two peakgroups.
#' @param stratify_regex NULL for no stratification or a string regular
#'   expression indicating categories which should drive peak splitting
#'   (e.g., batch or date).
#' @inheritParams aggregate_scans
#'
#' @returns a tibble containing two variables:
#' \itemize{
#'  \item{\code{groupId_old} - current groupIds in the mzroll_db_con}
#'  \item{\code{groupId_new} - updated groupIds in the mzroll_db_con}
#'  }
#'
#' @details
#' To be conservative, identified split peak pairs must satisfy all of the following conditions:
#' \itemize{
#'   \item{\code{mass agreement} - within the mass tolerance specified in the
#'     \code{clamr_config}}
#'   \item{\code{retention time agreement} - RTs are within
#'     \code{max_rt_deviation} of one another}
#'   \item{\code{mutual exclusivity} - log-abundances must be anti-correlated
#'     beyond \code{anticorrelation_co} and inter-quartile range of the signal
#'     split fracitons (fraction of signal in A versus B: A / (A + B)) above
#'     \code{signal_frac_IQR_co}.}
#'   \item{\code{batch-driven} [optional] - if batches are provided (using
#'     \code{stratify_regex}) then signal fraction variation should be explained
#'     by batches.}
#'   \item{\code{fragmentation agreement} - fragmentation spectra are correlated above \code{spectra_corr_co}}
#' }
#'
#' @export
identify_split_peaks <- function(mzroll_db_con, clamr_config, max_rt_deviation = 5, ic_floor = 2^10,
                                 anticorrelation_co = -0.1, signal_frac_IQR_co = 0.75, spectra_corr_co = 0.8,
                                 stratify_regex = NULL,
                                 n_top_spectra_summed = 3L, quality_weights = c("purity" = 2, "quality" = 1)) {
  debugr::dwatch(msg = "started identify_split_peaks. [clamr<processing.R>::identify_split_peaks.]\n")

  require_tolerances(clamr_config, c(1L, 2L))

  stopifnot(class(max_rt_deviation) == "numeric", length(max_rt_deviation) == 1, max_rt_deviation >= 0)
  stopifnot(class(anticorrelation_co) == "numeric", length(anticorrelation_co) == 1, anticorrelation_co <= 1, anticorrelation_co >= -1)
  stopifnot(class(ic_floor) == "numeric", length(ic_floor) == 1, ic_floor >= 0)
  stopifnot(any(c("NULL", "character") %in% class(stratify_regex)))
  if (class(stratify_regex) == "character") {
    stopifnot(length(stratify_regex) == 1)

    mzroll_samples <- tbl(mzroll_db_con, "samples") %>%
      dplyr::select(sampleId, name) %>%
      dplyr::collect() %>%
      dplyr::mutate(stratify_category = stringr::str_extract(name, stratify_regex))

    invalid_categories <- sum(is.na(mzroll_samples$stratify_category))
    if (invalid_categories != 0) {
      stop(invalid_categories, " were not matched by stratify_regex; either set stratify_regex to NULL or ensure that a pattern in all sample names is detected by stratify_regex")
    }

    n_statify_categories <- length(unique(mzroll_samples$stratify_category))
    if (n_statify_categories == 1) {
      warning("the provided stratify_regex did not distinguish samples into different categories")
    } else {
      message(n_statify_categories, " stratification categories found")
    }
  }

  # find pairs of peakgroups which are within mass-tolerance of one another
  peakset <- extract_peakset(mzroll_db_con)

  # compare all peakgroups by all peakgroups to find pairs of peakgroups within mass tolerance of one another

  variable_tolerances <- tibble::tribble(
    ~variable, ~tolerance, ~relative_or_absolute,
    "mean_peakMz", clamr_config$MS1tol$tol, clamr_config$MS1tol$absolute_or_relative
  )

  peakgroup_mass_matches <- join_matching_standards(
    peakset$peakgroups %>%
      dplyr::select(
        groupId_1 = groupId,
        mean_peakMz = mean_peakMz,
        mean_peakRt_1 = mean_peakRt
      ),
    peakset$peakgroups %>%
      dplyr::select(
        groupId_2 = groupId,
        mean_peakMz = mean_peakMz,
        mean_peakRt_2 = mean_peakRt
      ),
    variable_tolerances
  ) %>%
    # remove self-self and remove the redundant pair (keep group1:group2 and remove group2:group1)
    dplyr::filter(
      groupId_1 < groupId_2,
      # matches within a feasible rt deviation
      abs(mean_peakRt_1 - mean_peakRt_2) <= max_rt_deviation
    )

  candidate_tracking <- glue::glue("MZ consistent, RT within {max_rt_deviation} minutes - {nrow(peakgroup_mass_matches)}")

  # find the abundances of all pairs of feasible split peaks across all samples
  candidate_splits_abundances <- peakgroup_mass_matches %>%
    dplyr::select(groupId_1, groupId_2) %>%
    tidyr::crossing(peakset$peaks %>% dplyr::distinct(sampleId)) %>%
    dplyr::left_join(peakset$peaks %>%
      dplyr::select(groupId_1 = groupId, sampleId, peakAreaTop_1 = peakAreaTop),
    by = c("groupId_1", "sampleId")
    ) %>%
    dplyr::left_join(peakset$peaks %>%
      dplyr::select(groupId_2 = groupId, sampleId, peakAreaTop_2 = peakAreaTop),
    by = c("groupId_2", "sampleId")
    ) %>%
    # replace NAs with zeros
    dplyr::mutate(
      peakAreaTop_1 = ifelse(is.na(peakAreaTop_1), 0, peakAreaTop_1),
      peakAreaTop_2 = ifelse(is.na(peakAreaTop_2), 0, peakAreaTop_2),
      signal_frac_minor = dplyr::case_when(
        peakAreaTop_1 == 0 & peakAreaTop_2 == 0 ~ NA_real_,
        stats::median(peakAreaTop_1) > stats::median(peakAreaTop_2) ~ peakAreaTop_2 / (peakAreaTop_1 + peakAreaTop_2),
        TRUE ~ peakAreaTop_1 / (peakAreaTop_1 + peakAreaTop_2)
      )
    ) %>%
    # floor signals before moving to log-space
    dplyr::mutate(
      peakAreaTop_1 = ifelse(peakAreaTop_1 < ic_floor, ic_floor, peakAreaTop_1),
      peakAreaTop_2 = ifelse(peakAreaTop_2 < ic_floor, ic_floor, peakAreaTop_2),
      log2_ic_1 = log2(peakAreaTop_1),
      log2_ic_2 = log2(peakAreaTop_2)
    )

  # splits identified based on abundance relationships between the pair of peakgroups
  splits_summaries <- candidate_splits_abundances %>%
    dplyr::group_by(groupId_1, groupId_2) %>%
    dplyr::summarize(
      signal_frac_IQR = stats::IQR(signal_frac_minor, na.rm = TRUE),
      group_corr = stats::cor(log2_ic_1, log2_ic_2)
    ) %>%
    dplyr::ungroup() %>%
    # filter to
    dplyr::filter(
      signal_frac_IQR > signal_frac_IQR_co, # many near zero and near one values of signal frac results in an elevated IQR
      group_corr < anticorrelation_co
    ) # mutually exclusive splitting should induce a negative correlation

  candidate_tracking <- c(candidate_tracking, glue::glue("logIC < {anticorrelation_co}, signal_frac_IQR > {signal_frac_IQR_co} - {nrow(splits_summaries)}"))

  if (class(stratify_regex) == "character") {

    # identify cases where mutual exclusive signals are driven by batch effects
    # where batch effects are specified by the categories captured by the stratify_category

    candidate_splits_signal_stratification <- candidate_splits_abundances %>%
      dplyr::semi_join(splits_summaries, by = c("groupId_1", "groupId_2")) %>%
      dplyr::left_join(mzroll_samples, by = "sampleId") %>%
      tidyr::nest(group_data = -c("groupId_1", "groupId_2")) %>%
      dplyr::mutate(stratification_summary = furrr::future_map(group_data, stratify_signal_frac)) %>%
      dplyr::select(-group_data) %>%
      tidyr::unnest(stratification_summary) %>%
      dplyr::filter(is.na(p.value) | p.value < 0.01)

    candidate_tracking <- c(candidate_tracking, glue::glue("signal_frac variation driven by regex {stratify_regex} - {nrow(candidate_splits_signal_stratification)}"))

    splits_summaries <- splits_summaries %>%
      dplyr::semi_join(candidate_splits_signal_stratification, by = c("groupId_1", "groupId_2"))
  }

  # determine a consensus MS2 spectrum for each peakgroup
  consensus_peakgroup_spectra <- join_peaks_to_scans(mzroll_db_con) %>%
    dplyr::filter(!is.na(groupId)) %>%
    dplyr::filter(groupId %in% c(
      splits_summaries$groupId_1,
      splits_summaries$groupId_2
    )) %>%
    tidyr::nest(group_scans = -groupId) %>%
    dplyr::mutate(aggregated_scans = furrr::future_map(group_scans, aggregate_scans,
      MS2tol = clamr_config$MS2tol,
      n_top_spectra_summed = n_top_spectra_summed,
      quality_weights = quality_weights
    )) %>%
    dplyr::select(-group_scans) %>%
    tidyr::unnest_wider(aggregated_scans)

  # compare the fragmentation spectra of all candidate split pairs
  split_fragmentation_agreement <- splits_summaries %>%
    dplyr::left_join(consensus_peakgroup_spectra %>%
      dplyr::select(
        groupId_1 = groupId,
        fragmentationData_1 = fragmentationData_summary
      ),
    by = "groupId_1"
    ) %>%
    dplyr::left_join(consensus_peakgroup_spectra %>%
      dplyr::select(
        groupId_2 = groupId,
        fragmentationData_2 = fragmentationData_summary
      ),
    by = "groupId_2"
    ) %>%
    # consider pairs with fragmentation data for both peakgroups
    dplyr::filter(
      !purrr::map_lgl(fragmentationData_1, is.null),
      !purrr::map_lgl(fragmentationData_2, is.null)
    ) %>%
    # compare all scores
    dplyr::mutate(scores = purrr::map2(fragmentationData_1,
      fragmentationData_2,
      score_fragmentation_similarity,
      frag_similarity_methods = "robust_cosine",
      MS2tol = clamr_config$MS2tol
    )) %>%
    dplyr::select(-fragmentationData_1, -fragmentationData_2) %>%
    tidyr::unnest(scores) %>%
    dplyr::filter(method_score > spectra_corr_co)

  candidate_tracking <- c(candidate_tracking, glue::glue("MS2 robust correlation > {spectra_corr_co} - {nrow(split_fragmentation_agreement)}"))

  message(paste(candidate_tracking, collapse = "\n"))

  return(split_fragmentation_agreement %>%
    dplyr::select(groupId_1, groupId_2))
}

#' Find Unique Peakgroups
#'
#' For a set of peaks with similar m/z group peaks which have very similar rt
#' and correlated abundances across samples.
#'
#' @param mz_set_data a tibble containing groupIds
#' @param peaks a tibble containing all peaks in a mzrollDB
#' @inheritParams remove_redundant_peakgroups
#'
#' @return tibble of groupId and cluster integer
find_unique_peakgroups <- function(mz_set_data, peaks, rt_grouping_tol, group_corr_cutoff) {
  debugr::dwatch(msg = "Started find_unique_peakgroups. [clamr<processing.R>::find_unique_peakgroups]\n")

  stopifnot(class(rt_grouping_tol) == "numeric", length(rt_grouping_tol) == 1, rt_grouping_tol >= 0)
  stopifnot(class(group_corr_cutoff) == "numeric", length(group_corr_cutoff) == 1, group_corr_cutoff >= 0, group_corr_cutoff <= 1)

  reduced_peak_quants <- peaks %>%
    dplyr::semi_join(mz_set_data, by = "groupId") %>%
    dplyr::select(groupId, sampleId, ic = peakAreaTop)

  debugr::dwatch(msg = "determined reduced_peak_quants. [clamr<processing.R>::find_unique_peakgroups]\n")

  # If this is not handled as a special case, all groupIds will
  # be filtered out, which will produce a tibble with 0 rows.

  if (mz_set_data$groupId %>% length() == 1) {
    output <- tibble::tibble(groupId = mz_set_data$groupId, cluster = 1)
    return(output)
  }

  # find correlation between all groups
  groupId_correlations <- expand.grid(
    groupId_1 = mz_set_data$groupId, # TO DO: replace with dplyr::crossing
    groupId_2 = mz_set_data$groupId,
    stringsAsFactors = FALSE
  ) %>%
    # only look at upper-triangular entries
    dplyr::filter(groupId_1 < groupId_2) %>%
    dplyr::left_join(reduced_peak_quants %>%
      dplyr::rename(groupId_1 = groupId, ic_1 = ic),
    by = "groupId_1"
    ) %>%
    dplyr::left_join(reduced_peak_quants %>%
      dplyr::rename(groupId_2 = groupId, ic_2 = ic),
    by = c("sampleId", "groupId_2")
    ) %>%
    dplyr::group_by(groupId_1, groupId_2) %>%
    dplyr::summarize(corr = stats::cor(ic_1, ic_2, use = "pairwise.complete.obs")) %>%
    dplyr::ungroup()

  debugr::dwatch(msg = "determined groupId_correlations. [clamr<processing.R>::find_unique_peakgroups]\n")
  # print(groupId_correlations)

  # find pairs of peakgroups which agree based on correlation and rt
  groupId_matches <- groupId_correlations %>%
    dplyr::mutate(corr_match = corr > group_corr_cutoff) %>%
    dplyr::left_join(mz_set_data %>% dplyr::select(groupId_1 = groupId, rt_1 = mean_peakRt), by = "groupId_1") %>%
    dplyr::left_join(mz_set_data %>% dplyr::select(groupId_2 = groupId, rt_2 = mean_peakRt), by = "groupId_2") %>%
    dplyr::mutate(
      rt_match = abs(rt_1 - rt_2) < rt_grouping_tol,
      pair_match = corr_match & rt_match
    ) %>%
    dplyr::select(groupId_1, groupId_2, pair_match) %>%
    dplyr::filter(pair_match)

  debugr::dwatch(msg = "determined groupId_matches. [clamr<processing.R>::find_unique_peakgroups]\n")

  # use a greedy approach to group peakgroups into sparse subgraphs clusters
  group_clusters <- groupId_matches %>%
    igraph::graph_from_data_frame(directed = FALSE) %>%
    igraph::clusters() %>%
    igraph::membership() %>%
    {
      tibble::tibble(groupId = as.integer(names(.)), cluster = as.numeric(unname(.)))
    }

  debugr::dwatch(msg = "determined group_clusters. [clamr<processing.R>::find_unique_peakgroups]\n")

  max_defined_cluster <- ifelse(nrow(group_clusters) == 0, 0, max(group_clusters$cluster))

  # find singletons
  disjoint_groups <- mz_set_data %>%
    dplyr::anti_join(group_clusters, by = "groupId") %>%
    dplyr::select(groupId)

  debugr::dwatch(msg = "determined disjoint_groups. [clamr<processing.R>::find_unique_peakgroups]\n")

  if (nrow(disjoint_groups) == 0) {
    output <- group_clusters
  } else {
    disjoint_groups <- disjoint_groups %>%
      # singletons belong to own cluster
      dplyr::mutate(cluster = max_defined_cluster + 1:dplyr::n())

    output <- do.call(dplyr::bind_rows, list(group_clusters, disjoint_groups))
  }

  # print("find_unique_peakgroups output:")
  # print(output)

  output
}

#' Find M/Z Jumps
#'
#' Seperate an ordered vector of M/Z value into distinct groups with similar M/Z by finding gaps of greater than \code{cutoff}
#'
#' @param sorted_mzs numeric vector of m/z values
#' @param mass_accuracy_tolerance
#' \itemize{
#'   \item{If \code{absolute_or_relative} == "relative": minimum fractional mass tolerance (m/z * tol) between consecutive m/z values that indicates to start a new M/Z group}
#'   \item{If \code{absolute_or_relative} == "absolute": minimum mass tolerance (Da) between consecutive m/z values that indicates to start a new M/Z group}
#' }
#' @param absolute_or_relative Indicate whether the specified \code{mass_accuracy_tolerance} is "relative" or "absolute"
#' @param n_max_gap_inconsistencies Maximum number of observations of an mz-value to ignore when find breaks between common mzs (i.e., group A & B are seperate if there are <= n_max_gap_inconsistencies mzs within the tolerance between them)
#' @param collapse
#' \itemize{
#'   \item{If \code{collapse} == TRUE return mappings between each unique sorted_mz and mz_set}
#'   \item{If \code{collapse} == FALSE return a vector of mz_sets aligned to the sorted_mzs}
#' }
#'
#' @seealso \code{\link{find_density_minima}}
#'
#' @return a tibble containing paired mz and mz_set (mz groupings)
#'
#' @export
find_mz_jumps <- function(sorted_mzs, mass_accuracy_tolerance, absolute_or_relative, n_max_gap_inconsistencies = 0L, collapse = TRUE) {
  stopifnot(class(absolute_or_relative) == "character", length(absolute_or_relative) == 1, absolute_or_relative %in% c("relative", "absolute"))
  stopifnot(class(n_max_gap_inconsistencies) == "integer", length(n_max_gap_inconsistencies) == 1, n_max_gap_inconsistencies >= 0)
  stopifnot(class(collapse) == "logical", length(collapse) == 1, collapse %in% c(TRUE, FALSE))

  # identify jumps between blocks of m/z values to separate all peaks into clusters (which may contain multiple distinct species)

  mz_jumps <- rle(sorted_mzs) %>%
    {
      tibble::tibble(lengths = .$lengths, values = .$values)
    }

  if (n_max_gap_inconsistencies == 0) {
    mz_jumps$mzdiff_R <- c(NA, diff(mz_jumps$values))
  } else {
    # find sets of mzs which can be seperated allowing for skipping over a few observations

    # for each position, we want to determine how many indecies back we can go back (j) before we cross "n_max_gap_inconsistencies" observations
    # at this bound, we then test whether the differences between v[i] - v[i-j] > tolerance

    traceback_length <- rep(0, nrow(mz_jumps))
    traceback_cumsum <- rep(0, nrow(mz_jumps))
    for (i in 1:n_max_gap_inconsistencies) {
      traceback_cumsum[(i + 1):nrow(mz_jumps)] <- traceback_cumsum[(i + 1):nrow(mz_jumps)] + mz_jumps$lengths[1:(nrow(mz_jumps) - i)]
      traceback_cumsum[1:i] <- Inf
      traceback_length[traceback_cumsum <= n_max_gap_inconsistencies] <- i
    }
    mz_jumps$traceback_length <- traceback_length + 1
    # ensure that traceback stops at 1
    for (i in 1:(n_max_gap_inconsistencies + 1)) {
      mz_jumps$traceback_length[i] <- pmin(mz_jumps$traceback_length[i], i - 1)
    }
    mz_jumps$mzdiff_R <- mz_jumps$values - mz_jumps$values[c(1:nrow(mz_jumps)) - mz_jumps$traceback_length]

    # same process but look forward so that we aren't splitting groups that are similar to mzs which follow

    forwardsearch_length <- rep(0, nrow(mz_jumps))
    forwardsearch_cumsum <- rep(0, nrow(mz_jumps))
    for (i in 1:n_max_gap_inconsistencies) {
      forwardsearch_cumsum[1:(nrow(mz_jumps) - i)] <- forwardsearch_cumsum[1:(nrow(mz_jumps) - i)] + mz_jumps$lengths[(1 + i):nrow(mz_jumps)]
      forwardsearch_cumsum[nrow(mz_jumps):(nrow(mz_jumps) - i)] <- Inf
      forwardsearch_length[forwardsearch_cumsum <= n_max_gap_inconsistencies] <- i
    }
    mz_jumps$forwardsearch_length <- forwardsearch_length + 1
    # ensure that traceback stops at 1
    for (i in 1:(n_max_gap_inconsistencies + 1)) {
      mz_jumps$forwardsearch_length[nrow(mz_jumps) - i + 1] <- pmin(mz_jumps$forwardsearch_length[nrow(mz_jumps) - i + 1], i - 1)
    }
    mz_jumps$mzdiff_F <- mz_jumps$values[c(1:nrow(mz_jumps)) + mz_jumps$forwardsearch_length] - mz_jumps$values
  }

  # find unambiguously seperable sets of mzs
  if (absolute_or_relative == "relative") {
    mz_jumps$bin_start_R <- c(TRUE, ifelse(mz_jumps$mzdiff_R[-1] < mz_jumps$values[-1] * mass_accuracy_tolerance, FALSE, TRUE))
  } else {
    mz_jumps$bin_start_R <- c(TRUE, ifelse(mz_jumps$mzdiff_R[-1] < mass_accuracy_tolerance, FALSE, TRUE))
  }

  if (n_max_gap_inconsistencies != 0) {
    if (absolute_or_relative == "relative") {
      mz_jumps$bin_start_F <- c(TRUE, ifelse(mz_jumps$mzdiff_F[-1] < mz_jumps$values[-1] * mass_accuracy_tolerance, FALSE, TRUE))
    } else {
      mz_jumps$bin_start_F <- c(TRUE, ifelse(mz_jumps$mzdiff_F[-1] < mass_accuracy_tolerance, FALSE, TRUE))
    }
    mz_jumps$new_start <- c(TRUE, mz_jumps$bin_start_R[1:(nrow(mz_jumps) - 1)] == FALSE & mz_jumps$bin_start_R[2:nrow(mz_jumps)] == TRUE)
    mz_jumps$previous_start <- c(TRUE, mz_jumps$bin_start_R[1:(nrow(mz_jumps) - 1)] == TRUE & mz_jumps$bin_start_R[2:nrow(mz_jumps)] == TRUE & mz_jumps$bin_start_F[1:(nrow(mz_jumps) - 1)] == TRUE)
    mz_jumps$bin_start_consensus <- mz_jumps$new_start | mz_jumps$previous_start
    mz_jumps$mz_set <- cumsum(mz_jumps$bin_start_consensus)
  } else {
    mz_jumps$mz_set <- cumsum(mz_jumps$bin_start_R)
  }

  if (collapse == TRUE) {
    mz_jumps %>%
      dplyr::select(mz = values, mz_set)
  } else if (collapse == FALSE) {
    rep(mz_jumps$mz_set, times = mz_jumps$lengths)
  } else {
    stop('invalid "collapse" value')
  }
}

#' Returns a vector of grouped scans -- each group is defined by a separation of less than tolerances
#'
#' @description group sequential scans if scan[i] + tolerance[i] + tolerance[i + 1] > scan[i + 1]
#'
#' @param scans a vector of scan numbers
#' @param tolerances tolerances for grouping
#'
#' @return integer values of scan groups
find_scan_jumps <- function(scans, tolerances) {
  stopifnot(length(scans) == length(tolerances))

  if (length(scans) == 1) {
    return(1)
  }
  c(1, cumsum(scans[-length(scans)] + tolerances[-length(scans)] + tolerances[-1] < scans[-1]) + 1)
}

#' Smear Detection
#'
#' Identify and remove large background smears within a single precursor m/z.
#'
#' @param mz_sets tibble containing groups of peaks with a characteristic m/z.
#' @param rt_bin_cutoff_fraction within a given m/z window (mz_sets), filter the window if MS2 events are in more than this fraction of the 100 RT bins; a value of 1 filters nothing.
#' @param rescue_peaks TRUE/FALSE; if TRUE, within "smear" m/z regions, rescue some RT intervals which contain an IC pileup relative to adjacent RT bins; if FALSE, discard all data within a "smear" m/z region.
#' @param create_plots TRUE/FALSE; create a plot which shows the number of RT bins and total MS2 events along with \code{rt_bin_cutoff_fraction} as a "smear" cutoff
#'
#' @return tibble which is the same as the input only with smear m/zs removed.
#'
#' @export
smear_detection <- function(mz_sets, rt_bin_cutoff_fraction = 0.5, rescue_peaks = TRUE, create_plots = TRUE) {
  stopifnot(rt_bin_cutoff_fraction > 0, rt_bin_cutoff_fraction <= 1)
  if (rt_bin_cutoff_fraction == 1) {
    print("No smear filtering")
    return(mz_sets)
  }

  stopifnot(class(rescue_peaks) == "logical", rescue_peaks %in% c(TRUE, FALSE))
  stopifnot(class(create_plots) == "logical", create_plots %in% c(TRUE, FALSE))

  rt_bins <- (seq(min(mz_sets$rt), max(mz_sets$rt), length.out = 101)[1:100] + seq(min(mz_sets$rt), max(mz_sets$rt), length.out = 101)[2:101]) / 2

  mz_sets$rt_bin <- sapply(mz_sets$rt, function(an_rt) {
    which.min(abs(rt_bins - an_rt))
  })

  n_samples <- length(unique(mz_sets$sampleId))

  # determine the number of peaks (and their IC) within each RT bin

  mz_set_rt_density <- mz_sets %>%
    dplyr::group_by(mz_set, precursor_charge, rt_bin) %>%
    dplyr::summarize(
      n = dplyr::n(),
      ICsum = sum(as.numeric(precursor_ic))
    )

  # separate smear and non-smear

  mz_set_bin_summary <- mz_set_rt_density %>%
    dplyr::group_by(mz_set, precursor_charge) %>%
    dplyr::summarize(
      n_ms2_total = sum(n),
      n_bins = dplyr::n()
    )

  if (create_plots) {
    print(
      ggplot(
        mz_set_bin_summary %>%
          dplyr::inner_join(mz_set_bin_summary %>%
            dplyr::ungroup() %>%
            dplyr::count(precursor_charge) %>%
            dplyr::arrange(desc(n)) %>%
            dplyr::slice(1:4),
          by = "precursor_charge"
          ),
        aes(x = n_ms2_total, y = n_bins)
      ) +
        geom_point() +
        facet_wrap(~precursor_charge) +
        ggtitle("MS2 spectra for a given precursor mass") +
        # scale_x_log10("Number of MS2 spectra") +
        scale_y_continuous("Retention times with 1+ MS2 spectra (of 100 bins)") +
        theme_minimal()
    )
  }

  # call smear if more than rt_bin_cutoff_fraction fraction of retention time bins contain an MS2

  smeared_precursors <- mz_set_bin_summary %>%
    dplyr::mutate(is_smear = ifelse(n_bins <= 100 * rt_bin_cutoff_fraction, FALSE, TRUE))

  if (rescue_peaks == FALSE) {
    # early out if everythin within smear is discarded
    return(
      mz_sets %>%
        # clean background
        dplyr::semi_join(smeared_precursors %>%
          dplyr::filter(!is_smear),
        by = c("mz_set", "precursor_charge")
        )
    )
  }

  # identify peaks within smears

  mz_set_rt_density_smears <- mz_set_rt_density %>%
    dplyr::semi_join(smeared_precursors %>%
      dplyr::filter(is_smear),
    by = c("mz_set", "precursor_charge")
    )

  rt_bin_IC_enrichments <- mz_set_rt_density_smears %>%
    # find the ion count of all precursors -2 : +2 RT bins from RT_bin
    # -1 : +1 will be combined with RT_bin's signal in case a compound spans a bin's bound
    # -2 & +2 will be used to capture the peaks background (since it is a smear afterall)
    dplyr::left_join(mz_set_rt_density_smears %>%
      dplyr::mutate(rt_bin = rt_bin + 2) %>%
      dplyr::select(mz_set, precursor_charge, rt_bin, `n-2 IC` = ICsum),
    by = c("mz_set", "precursor_charge", "rt_bin")
    ) %>%
    dplyr::left_join(mz_set_rt_density_smears %>%
      dplyr::mutate(rt_bin = rt_bin + 1) %>%
      dplyr::select(mz_set, precursor_charge, rt_bin, `n-1 IC` = ICsum),
    by = c("mz_set", "precursor_charge", "rt_bin")
    ) %>%
    dplyr::left_join(mz_set_rt_density_smears %>%
      dplyr::mutate(rt_bin = rt_bin - 1) %>%
      dplyr::select(mz_set, precursor_charge, rt_bin, `n+1 IC` = ICsum),
    by = c("mz_set", "precursor_charge", "rt_bin")
    ) %>%
    dplyr::left_join(mz_set_rt_density_smears %>%
      dplyr::mutate(rt_bin = rt_bin - 2) %>%
      dplyr::select(mz_set, precursor_charge, rt_bin, `n+2 IC` = ICsum),
    by = c("mz_set", "precursor_charge", "rt_bin")
    ) %>%
    dplyr::rowwise() %>%
    # set a minimum ion count (conservative and avoids NaNs)
    dplyr::mutate_each(dplyr::funs(sum(c(., 300), na.rm = TRUE)), ICsum, `n-2 IC`, `n-1 IC`, `n+1 IC`, `n+2 IC`) %>%
    # compare signal under query bin and adjacent bin to 2 bins away
    dplyr::mutate(
      `+- 1 bin` = sum(`n-1 IC`, `n+1 IC`),
      `+- 2 bin` = sum(`n-2 IC`, `n+2 IC`) * 3 / 2,
      bin_IC_enrichment = (ICsum + `+- 1 bin`) / `+- 2 bin`
    ) %>%
    dplyr::ungroup()

  rt_bin_peaks <- rt_bin_IC_enrichments %>%
    dplyr::left_join(mz_sets %>%
      dplyr::group_by(mz_set, precursor_charge) %>%
      dplyr::summarize(ms1_mz = stats::median(precursor_mz)),
    by = c("mz_set", "precursor_charge")
    ) %>%
    # filter to identify peaks based on a 5-fold enrichment of signal relative to background
    # number of MS2s doesn't seem to track with peak quality
    dplyr::filter(bin_IC_enrichment > 5, n > 1)

  if (nrow(rt_bin_peaks) != 0) {
    smear_peaks <- lapply(1:nrow(rt_bin_peaks), function(i) {
      rt_indices <- (rt_bin_peaks$rt_bin[i] - 2):(rt_bin_peaks$rt_bin[i] + 2)
      rt_indices <- rt_indices[rt_indices > 0]
      data.frame(mz_set = rt_bin_peaks$mz_set[i], precursor_charge = rt_bin_peaks$precursor_charge[i], rt_bin = rt_indices)
    }) %>%
      dplyr::bind_rows()

    mz_sets_smear_peaks <- mz_sets %>%
      dplyr::semi_join(smear_peaks, by = c("mz_set", "precursor_charge", "rt_bin"))
  } else {
    mz_sets_smear_peaks <- NULL
  }

  mz_sets %>%
    # clean background (non-smear)
    dplyr::semi_join(smeared_precursors %>%
      dplyr::filter(!is_smear),
    by = c("mz_set", "precursor_charge")
    ) %>%
    # peaks within smears
    dplyr::bind_rows(mz_sets_smear_peaks) %>%
    dplyr::select(-rt_bin)
}

#' Reextracted Peak Intensity Set
#'
#' Produce a table of reextracted peak intensities based on a short
#' list of curated peak groups
#'
#' @param curated_mzrolldb
#'     mzrollDB file containing manually verified peakgroups.
#'     peaks and peakgroups from this table are used in this method.
#'
#' @param full_mzrolldb
#'     mzrollDB file containing original rt_update_key RT alignment table.
#'     This table is used for transforming sample_rt <--> aligned_rt
#'
#' @param sample_folder
#'     folder containing all raw samples to apply reextraction.
#'
#' @param min_mz_delta
#'     adjust min and max m/z so that the difference is at least this value.
#'
#' @param min_rt_delta
#'     adjust min and max rt so that the difference is at least this value.
#'
#' @param verbose
#'     flag for explicit printing
#'
#' @param debug
#'     flag for debugging
#'
#' @return intensities_table
#'
#'    Each row indicates quant values for a single peak (an annotated feature in a single sample.)
#'    Both aligned RT coordinates and sample RT coordinates are given.
#'
#'    \item{sample: }{sample name}
#'    \item{sampleId: }{ID number for sample in mzrolldb files}
#'    \item{groupId: }{group ID number from curated_mzrolldb file}
#'    \item{compoundName: }{String representation of annotated feature}
#'    \item{intensity: }{curated_mzrollDB-recorded intensity (peakAreaTop), or reextracted intensity}
#'    \item{is_identified: }{if TRUE, using quant from curated_mzrolldb. Else, reextract from sample file}
#'    \item{mzmin: }{minimum m/z associated with peak}
#'    \item{mzmax: }{maximum m/z associated with peak}
#'    \item{aligned_rtmin: }{minimum RT value associated with peak, in aligned RT space. Relies on rt_update_key table in full_mzrolldb}
#'    \item{aligned_rtmax: }{maximum RT value associated with peak, in aligned RT space. Relies on rt_update_key table in full_mzrolldb}
#'    \item{sample_rtmin: }{minimum RT value associated with peak, in mzML/mzXML sample space.}
#'    \item{sample_rtmax: }{maximum RT value associated with peak, in mzML/mzXML sample space.}
#' @export
reextracted_compound_intensities <- function(curated_mzrolldb,
                                             full_mzrolldb,
                                             sample_folder,
                                             min_mz_delta = 0.01,
                                             min_rt_delta = 0.1,
                                             verbose = TRUE,
                                             debug = FALSE) {
  if (!rlang::is_installed(("mzkitcpp"))) {
    stop("Error: mzkitcpp must be installed to use this function.")
  }

  # Need RT alignment information for all samples for accurate reextraction
  full_mzroll_db_con <- open_mzrolldb(full_mzrolldb)

  rt_alignment <- DBI::dbGetQuery(
    full_mzroll_db_con,
    "SELECT samples.name, rt_update_key.sampleId, rt_update_key.rt, rt_update_key.rt_update FROM rt_update_key
           INNER JOIN samples on rt_update_key.sampleId = samples.sampleId"
  )

  DBI::dbDisconnect(full_mzroll_db_con)

  # Other aspects of the search come from the manually curated mzrollDB
  mzroll_db_con <- open_mzrolldb(curated_mzrolldb)

  peaks_data <- DBI::dbGetQuery(
    mzroll_db_con,
    "SELECT peaks.peakId,peaks.groupId, peaks.sampleId, samples.name, peaks.rtmin, peaks.rtmax, peaks.mzmin, peaks.mzmax,
            peaks.peakArea, peaks.peakAreaCorrected, peaks.peakAreaTop, peaks.peakIntensity FROM peaks
           INNER JOIN samples on peaks.sampleId = samples.sampleId;"
  )

  peakgroups_data <- DBI::dbGetQuery(
    mzroll_db_con,
    "SELECT peakgroups.groupId, peakgroups.parentGroupId, peakgroups.compoundName, peakgroups.adductName, peakgroups.tagString FROM peakgroups
        WHERE
           peakgroups.searchTableName IS 'rumsDB' AND peakgroups.label LIKE '%g%'
        OR peakgroups.searchTableName IS 'clamDB' AND peakgroups.label LIKE '%g%'
        OR peakgroups.searchTableName IS 'Bookmarks';"
  )

  DBI::dbDisconnect(mzroll_db_con)

  name_to_sample_id <- rt_alignment %>%
    dplyr::select(name, sampleId) %>%
    unique()

  # Convert RT values to mzML/mzXML space prior to reextraction
  #
  # RT values for peaks and peakgroups corresponds to post-alignment, rt_update_key table "rt_update" column
  #
  # (Aligner::doSegmentedAligment() adjusts all scan RT values,
  # then peakdetector uses the updated scan RTs to form peaks and peakgroups)
  #
  # RT values in mzML / mzXML files correspond to pre-alignment, rt_update_key table "rt" column (see ProjectDB::saveAlignment())
  #
  # RT alignment is always carried out using the column named "rt" as the starting point, and projects onto "rt_update" space
  # So to invert RT transformation, need to switch column names.

  rt_alignment_inverted <- rt_alignment %>%
    dplyr::rename(rt_update = rt, rt = rt_update)

  peaks_rt_min_aligned <- mzkitcpp::update_rts(rt_alignment_inverted, peaks_data$rtmin, peaks_data$sampleId)
  peaks_rt_max_aligned <- mzkitcpp::update_rts(rt_alignment_inverted, peaks_data$rtmax, peaks_data$sampleId)

  peaks_data_w_sample_rt <- cbind(
    peaks_data,
    tibble::tibble(
      "sample_rtmin" = peaks_rt_min_aligned$updated_rts,
      "sample_rtmax" = peaks_rt_max_aligned$updated_rts
    )
  )

  # compound_list RTs are in aligned RT space
  compound_list <- dplyr::inner_join(peakgroups_data, peaks_data, by = c("groupId")) %>%
    dplyr::mutate(mz_delta = mzmax - mzmin) %>%
    dplyr::mutate(rt_delta = rtmax - rtmin) %>%
    dplyr::mutate(mz_median = mzmin + 0.5 * mz_delta) %>%
    dplyr::mutate(rt_median = rtmin + 0.5 * rt_delta) %>%
    # Adjust m/z and Rt windows to be at least as wide as the min_mz_delta / min_rt_delta
    dplyr::mutate(mzmin_adj = ifelse(mz_delta < min_mz_delta, mz_median - 0.5 * min_mz_delta, mzmin)) %>%
    dplyr::mutate(mzmax_adj = ifelse(mz_delta < min_mz_delta, mz_median + 0.5 * min_mz_delta, mzmax)) %>%
    dplyr::mutate(rtmin_adj = ifelse(rt_delta < min_rt_delta, rt_median - 0.5 * min_rt_delta, rtmin)) %>%
    dplyr::mutate(rtmax_adj = ifelse(rt_delta < min_rt_delta, rt_median + 0.5 * min_rt_delta, rtmax)) %>%
    # Mutate compound name to include isotope information, or unnamed parts, if applicable
    dplyr::mutate(
      compoundName =
        paste0(
          ifelse(!is.na(compoundName), compoundName, paste0(round(mz_median, digits = 2), "@", round(rt_median, digits = 2))),
          ifelse(tagString == "", "", paste0(" ", tagString)),
          ifelse(is.na(adductName) | is.na(compoundName), "", paste0(" ", adductName))
        )
    ) %>%
    # select columns of interest
    dplyr::select(groupId, compoundName, rtmin_adj, rtmax_adj, mzmin_adj, mzmax_adj) %>%
    dplyr::rename(mzmin = mzmin_adj, mzmax = mzmax_adj, rtmin = rtmin_adj, rtmax = rtmax_adj) %>%
    # generate lenient bounds based on identified samples
    dplyr::group_by(groupId) %>%
    dplyr::mutate(mzmin = min(mzmin), mzmax = max(mzmax), rtmin = min(rtmin), rtmax = max(rtmax)) %>%
    dplyr::ungroup() %>%
    unique() %>%
    dplyr::rename(aligned_rtmin = rtmin, aligned_rtmax = rtmax)

  if (verbose) {
    cat(paste0("Retrieved ", nrow(compound_list), " annotated peak groups from curated mzroll db file \"", basename(curated_mzrolldb), "\"\n"))
  }

  sample_files <- list.files(sample_folder, pattern = "mzX?ML", full.names = TRUE)

  if (verbose) {
    cat(paste0("Discovered ", length(sample_files), " mzML/mzXML sample files for reextraction.\n"))
  }

  reextracted_compound_list <- tibble::tibble(
    "sample" = character(0),
    "sampleId" = integer(0),
    "groupId" = integer(0),
    "compoundName" = character(0),
    "intensity" = numeric(0),
    "is_identified" = logical(0),
    "mzmin" = numeric(0),
    "mzmax" = numeric(0),
    "aligned_rtmin" = numeric(0),
    "aligned_rtmax" = numeric(0),
    "sample_rtmin" = numeric(0),
    "sample_rtmax" = numeric(0)
  )

  compound_list_compounds <- compound_list %>% dplyr::select(groupId, compoundName)

  for (i in 1:length(sample_files)) {

    # Determine all sample-specific information
    sample_file <- sample_files[i]
    sample_name <- basename(sample_file)
    sample_id_tbl <- name_to_sample_id %>% dplyr::filter(name == sample_name)
    sample_id <- sample_id_tbl$sampleId[1]

    # Recover quant for samples that have quant values already
    identified_peaks <- peaks_data_w_sample_rt %>%
      dplyr::filter(name == sample_name) %>%
      dplyr::inner_join(compound_list_compounds, by = c("groupId")) %>%
      dplyr::rename(
        intensity = peakAreaTop,
        aligned_rtmin = rtmin,
        aligned_rtmax = rtmax,
        sample = name
      ) %>%
      dplyr::mutate(is_identified = TRUE) %>%
      dplyr::select(sample, sampleId, groupId, compoundName, intensity, is_identified, mzmin, mzmax, aligned_rtmin, aligned_rtmax, sample_rtmin, sample_rtmax)

    # Remove identified peaks from the list of compounds to search
    compound_list_ith_sample <- compound_list %>%
      dplyr::filter(!groupId %in% identified_peaks$groupId)

    if (nrow(identified_peaks) > 0) {
      compound_list_ith_sample <- compound_list %>%
        dplyr::filter(!compoundName %in% identified_peaks$compoundName)
    } else {
      compound_list_ith_sample <- compound_list
    }

    # add sample RTs to allow for reextraction
    compound_list_ith_sample_rtmin <- mzkitcpp::update_rts(rt_alignment_inverted, compound_list_ith_sample$aligned_rtmin, rep(sample_id, nrow(compound_list_ith_sample)))
    compound_list_ith_sample_rtmax <- mzkitcpp::update_rts(rt_alignment_inverted, compound_list_ith_sample$aligned_rtmax, rep(sample_id, nrow(compound_list_ith_sample)))

    compound_list_ith_sample_reextract <- cbind(
      compound_list_ith_sample,
      "rtmin" = compound_list_ith_sample_rtmin$updated_rts,
      "rtmax" = compound_list_ith_sample_rtmax$updated_rts
    )

    # Re-extract any unidentified compounds.
    # Reextracted values returned here are in mzXML/mzML RT space
    reextracted_peaks <- mzkitcpp::reextract_peaks(
      compound_list_ith_sample_reextract,
      sample_file,
      verbose = FALSE,
      debug = debug
    ) %>%
      dplyr::rename(
        sample_rtmin = rtmin,
        sample_rtmax = rtmax,
        compoundName = compound
      ) %>%
      dplyr::inner_join(name_to_sample_id, by = c("sample" = "name")) %>%
      dplyr::inner_join(compound_list_ith_sample_reextract, by = c("groupId")) %>%
      dplyr::rename(
        mzmin = mzmin.x,
        mzmax = mzmax.x,
        compoundName = compoundName.x
      ) %>%
      dplyr::mutate(is_identified = FALSE) %>%
      dplyr::select(sample, sampleId, groupId, compoundName, intensity, is_identified, mzmin, mzmax, aligned_rtmin, aligned_rtmax, sample_rtmin, sample_rtmax)

    ith_sample_intensities <- rbind(identified_peaks, reextracted_peaks)

    if (verbose) {
      cat(paste0(
        "Completed peak reextraction for sample ",
        i,
        "/",
        length(sample_files),
        ": \"",
        basename(sample_files[i]),
        "\""
      ))
    }
    reextracted_compound_list <- rbind(
      reextracted_compound_list,
      ith_sample_intensities
    )
  }

  reextracted_compound_list
}

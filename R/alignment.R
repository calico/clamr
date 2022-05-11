#' MS2 Driven Aligner
#'
#' @description Across a dataset, use MS1 and MS2 data to generate update function which can be used to transform retention times and mass accuracy in order to improve consistency.
#'
#' @inheritParams test_mzroll_db_con_schema
#' @inheritParams test_clamr_config
#' @param group_by_charge Require that peaks match in precursor charge in order to be a possible match -- discard all peaks with unknown charge.
#' @param peak_quality_cutoff Minimum quality for a scan's matching peak for the ion to be used for alignment.
#' @param maximum_mz_set_size The maximum number of possibly matching ms2 events that will be considered.
#' @param return_plot_data return data which can be used to plot alignment summaries with \code{\link{ms2_driven_aligner_plotting}}
#' @inheritParams find_consistent_ms2_fingerprints
#' @inheritParams estimate_rt_drift
#' @inheritParams build_clamr_config
#'
#' @details This function uses MS1 and MS2 similarity to group MS2 events from different samples together and then determines the extent to which
#' individual samples systematically deviate in either retention time or mass accuracy.
#'
#' @return an update mzroll_db_con where all features have had a function applied of the form RT_updated = RT_original + g(RT_original), where the aim of the g(RT_original) is to estimate retention time deviations between individual samples and a consensus sample.
#'
#' @examples
#' \dontrun{
#' # TO DO - add a smaller dataset so this would run fast enough for testing
#' library(dplyr)
#'
#' mzroll_db_con <- clamshell::get_clamr_assets("mzkit_peakdetector.mzrollDB") %>%
#'   mzroll_db_sqlite()
#'
#' clamr_config <- build_clamr_config(list(MS1tol = "10ppm", MS2tol = "20ppm"))
#' ms2_driven_aligner(mzroll_db_con, clamr_config)
#' }
#'
#' @export
ms2_driven_aligner <- function(mzroll_db_con,
                               clamr_config,
                               group_by_charge = FALSE,
                               peak_quality_cutoff = 0.5,
                               maximum_mz_set_size = 1000,
                               cosine_cutoff = 0.95,
                               sd_rt_resid_cutoff = 0.1,
                               spline_ridge_penalty = 200,
                               spline_degree = 4L,
                               return_plot_data = FALSE,
                               quietly = FALSE) {

  # confirm that data contains MS2 spectra and appropriate experimental parameters

  test_mzroll_db_con_schema(mzroll_db_con)
  require_tolerances(clamr_config, required_msLevels = c(1L, 2L))

  # test validity of other parameters
  stopifnot(class(group_by_charge) == "logical", length(group_by_charge) == 1, group_by_charge %in% c(TRUE, FALSE))
  stopifnot(class(peak_quality_cutoff) %in% c("logical", "numeric"), length(peak_quality_cutoff) == 1)
  if (peak_quality_cutoff == "logical") {
    stopifnot(is.na(peak_quality_cutoff))
  } else if (peak_quality_cutoff == "numeric") {
    stopifnot(peak_quality_cutoff >= 0, peak_quality_cutoff <= 1)
  }
  stopifnot(class(maximum_mz_set_size) %in% c("numeric", "integer"), length(maximum_mz_set_size) == 1, maximum_mz_set_size > 0)
  stopifnot(length(cosine_cutoff) == 1, cosine_cutoff >= 0, cosine_cutoff < 1)
  stopifnot(length(sd_rt_resid_cutoff) == 1, sd_rt_resid_cutoff > 0)
  stopifnot(length(spline_ridge_penalty) == 1, class(spline_ridge_penalty) == "numeric", spline_ridge_penalty >= 0)
  stopifnot(length(spline_degree) == 1, class(spline_degree) == "integer", spline_degree >= 3L)
  stopifnot(class(return_plot_data) == "logical", length(return_plot_data) == 1, return_plot_data %in% c(TRUE, FALSE))

  # Group peaks (drawn from all samples based on MS1 m/z similarity)
  # separate all peak groups to find large breaks in parent ion m/zs
  # this is fast as is but doesn't currently take advantage of both lists being ordered


  # join peaks to scans when there is a scan range and mz match for a given sample
  debugr::dwatch(msg = "Joining mzrolldb peaks table to mzrolldb scans table. [clamr<alignment.R>::ms2_driven_aligner]")

  scan_peak_join <- join_peaks_to_scans(mzroll_db_con)

  if (length(unique(scan_peak_join$sampleId)) <= 1) {
    stop("<2 sampleId detected: multiple samples are required for alignment")
  }

  # filter by quality score (optionally)

  if (!(is.na(peak_quality_cutoff))) {
    # only take scans which have been matched to a peak, and the peak quality is over the specified threshold
    filtered_scan_peak_join <- scan_peak_join %>%
      dplyr::filter(!is.na(quality)) %>%
      dplyr::filter(quality > peak_quality_cutoff)
  } else {
    filtered_scan_peak_join <- scan_peak_join
  }

  debugr::dwatch(msg = "Overwriting scan mz and rt with peak attributes (if available). [clamr<alignment.R>::ms2_driven_aligner]")
  # overwrite scan mz and rt with peak attributes if available

  filtered_scan_peak_join <- filtered_scan_peak_join %>%
    dplyr::mutate(
      rt = ifelse(is.na(peakRt), rt, peakRt),
      precursorMz = ifelse(is.na(peakMz), precursorMz, peakMz),
      precursorIc = ifelse(is.na(peakAreaTop), precursorIc, peakAreaTop)
    ) %>%
    # remove irrelevant variables
    dplyr::select(sampleId, scan, precursorCharge, rt, precursorMz, precursorIc, data)

  # either use precursor_charge for grouping or ignore precursor charge and allow for ions with unknown charge
  debugr::dwatch(msg = "Adjust scans based on precursor charge information (if available). [clamr<alignment.R>::ms2_driven_aligner]")

  if (group_by_charge) {
    unknown_charge_fraction <- sum(filtered_scan_peak_join$precursorCharge == 0) / nrow(filtered_scan_peak_join)
    if (unknown_charge_fraction > 0.05) {
      warning(round(unknown_charge_fraction * 100, 1), '% of MS2 spectra have an unknown precursor charge value
              and will not be used. Use "group_by_charge = FALSE" to included these peaks')
    }
    filtered_scan_peak_join <- filtered_scan_peak_join %>%
      dplyr::filter(precursorCharge != 0)
  } else {
    # zero out charges if they will not be used
    filtered_scan_peak_join <- filtered_scan_peak_join %>%
      dplyr::mutate(precursorCharge = -1L)
  }

  debugr::dwatch(msg = "Creating mz_sets. [clamr<alignment.R>::ms2_driven_aligner]")
  mz_sets <- filtered_scan_peak_join %>%
    dplyr::group_by(precursorCharge) %>%
    # within each charge group sort by precursor m/z
    dplyr::arrange(precursorMz) %>%
    # find and label groups with similar m/z
    dplyr::mutate(mz_set = find_mz_jumps(.$precursorMz,
      clamr_config$MS1tol$tol,
      clamr_config$MS1tol$absolute_or_relative,
      collapse = FALSE
    )) %>%
    dplyr::ungroup()

  # organize scans with matching precursorMz and precursorCharge into clusters with similar MS2 fragmentation profiles based on cosine similarity.
  debugr::dwatch(msg = "Organizing scans with matching precursorMz and precursorCharge into clusters with similar MS2 fragmentation profiles. [clamr<alignment.R>::ms2_driven_aligner]")

  all_ms2_groups <- mz_sets %>%
    tidyr::nest(mz_set_data = c(-mz_set, -precursorCharge)) %>%
    # filter groups that are greater than a maximum groupsize where calculating an O(n^2) distance matrix would be computationally infeasible
    dplyr::mutate(n_entries = purrr::map_int(mz_set_data, nrow)) %>%
    dplyr::filter(n_entries > 1 & n_entries < maximum_mz_set_size) %>%
    # for peaks with the same precursor (MS1) mass, group peaks based on MS2 similarity
    dplyr::mutate(reduced_groups = furrr::future_map(mz_set_data,
      group_similar_MS2,
      clamr_config = clamr_config,
      cosine_cutoff = cosine_cutoff
    )) %>%
    dplyr::filter(!purrr::map_lgl(reduced_groups, is.null)) %>%
    tidyr::unnest_wider(reduced_groups)

  # Summarize features to determine how reproducibly they are measured across all samples
  debugr::dwatch(msg = "Summarizing features to determine cross-sample reproducibility. [clamr<alignment.R>::ms2_driven_aligner]")

  ms2_groups <- all_ms2_groups %>%
    dplyr::select(-mz_set_data, -ms2_footprints) %>%
    tidyr::unnest(ms2_groups)

  ms2_peak_consistency <- ms2_groups %>%
    dplyr::distinct(sampleId, mz_set, precursorCharge, ms2_group) %>%
    dplyr::count(mz_set, precursorCharge, ms2_group) %>%
    dplyr::filter(n > 1) %>%
    # create a unique label for each group
    dplyr::mutate(compound_id = glue::glue("{mz_set}_{precursorCharge}_{ms2_group}"))

  # Groups of ostensibly matched compounds (compound_id) along with the sample they are from and retention time observed
  debugr::dwatch(msg = "Group matched compounds, associate with correspondign sample and observed RT. [clamr<alignment.R>::ms2_driven_aligner]")

  sample_matched_compounds <- ms2_groups %>%
    # remove compound groups filtered in ms2_peak_consistency
    dplyr::inner_join(ms2_peak_consistency, by = c("mz_set", "precursorCharge", "ms2_group")) %>%
    dplyr::select(sampleId, compound_id, precursorMz, rt)

  # Fitting of retention time deviation -- iterative estimation of "true" compound RT and sample-specific deviation between true and observed.
  debugr::dwatch(msg = "Fitting RT deviation. [clamr<alignment.R>::ms2_driven_aligner]")

  # Check existing update key, inverse transformation, if necessary
  original_rt_update_key_swapped <- dplyr::tbl(mzroll_db_con, dplyr::sql("SELECT sampleId, rt, rt_update FROM rt_update_key")) %>%
    dplyr::collect() %>%
    dplyr::rename(rt_update = rt, rt = rt_update)

  debugr::dwatch(msg = "Recovered original_rt_update_key_swapped. [clamr<alignment.R>::ms2_driven_aligner]")

  fit_RT_update <- estimate_rt_drift(sample_matched_compounds,
    previous_rt_update_key = original_rt_update_key_swapped,
    sd_rt_resid_cutoff = sd_rt_resid_cutoff,
    spline_ridge_penalty = spline_ridge_penalty,
    spline_degree = spline_degree,
    return_plot_data = return_plot_data
  )

  output <- list()
  output$rt_drift <- fit_RT_update$sample_gam_fits
  output$rt_sd <- fit_RT_update$rt_sd
  output$rt_update_key <- fit_RT_update$rt_update_key

  # Fitting ppm drift in mass accuracy
  debugr::dwatch(msg = "Fitting ppm drift in mass accuracy. [clamr<alignment.R>::ms2_driven_aligner]")

  if (clamr_config$MS1tol$absolute_or_relative == "relative") {
    fit_MZ_update <- estimate_mz_ppm_drift(sample_matched_compounds, return_plot_data = return_plot_data)

    output$MS1_drift <- fit_MZ_update$MS1_drift
    output$MS1_sd <- fit_MZ_update$MS1_sd
  } else {
    warning('"clamr_config$MS1tol$absolute_or_relative" can only be specified as relative (currently)
            when correcting for drift in mass accuracy\nno mass accuracy update information will be returned')
  }

  if (return_plot_data) {
    debugr::dwatch(msg = "Generating plot data to be returned. [clamr<alignment.R>::ms2_driven_aligner]")

    output$plot_data <- generate_rt_plot_data(
      mzroll_db_con,
      clamr_config,
      ms2_groups,
      ms2_peak_consistency,
      fit_RT_update,
      fit_MZ_update
    )
  }

  debugr::dwatch(msg = "Returning output. [clamr<alignment.R>::ms2_driven_aligner]")
  output
}

#' Estimate RT Drift
#'
#' @description Generate a mapping of observed sample retention times onto imputed (consistent) retention time.
#'
#' @details Alternate between estimating compound-specific retention time and sample-specific retention time deviations.
#'
#' @param sample_matched_compounds tibble containing at least sampleId, compound_id and rt.
#' @param previous_rt_update_key previous mapping of original RTs in raw files to RTs used in peaks and scans.
#' @param sd_rt_resid_cutoff Cutoff for excluding a compound based on the standard deviation of its residuals (rt - fitted rt): sd(resid)/range(rt) < \code{sd_rt_resid_cutoff}.
#' @param spline_ridge_penalty Ridge penalty on sample-level splines which model deviations between observed retention times and the consensus rt of a compound. Used for the H parameter in \link[mgcv]{gam}.
#' @param spline_degree Degree of spline used to estimate drift (integer).
#' @inheritParams ms2_driven_aligner
#'
#' @return a list containing sample specific gam fits and a sample independent fit which should be used to correct sample-level overfitting
estimate_rt_drift <- function(sample_matched_compounds, previous_rt_update_key = NULL, sd_rt_resid_cutoff = 0.1, spline_ridge_penalty = 200, spline_degree = 4L, return_plot_data = FALSE) {
  if (debugr::debugr_isActive()) {
    input_variables_msg <- paste0(
      "[clamr<alignment.R>::estimate_rt_drift] INPUT PARAMETERS:\n",
      "sample_matched_compounds= is.null? ", is.null(sample_matched_compounds),
      "\nprevious_rt_update_key= is.null? ", is.null(previous_rt_update_key),
      "\nsd_rt_resid_cutoff=", sd_rt_resid_cutoff,
      "\nspline_ridge_penalty=", spline_ridge_penalty,
      "\nspline_degree=", spline_degree,
      "\nreturn_plot_data=", return_plot_data
    )
    debugr::dwatch(msg = input_variables_msg)
  }

  stopifnot(length(sd_rt_resid_cutoff) == 1, class(sd_rt_resid_cutoff) %in% c("numeric", "integer"), sd_rt_resid_cutoff > 0)
  stopifnot(length(spline_ridge_penalty) == 1, class(spline_ridge_penalty) %in% c("numeric", "integer"), spline_ridge_penalty >= 0)
  stopifnot(length(spline_degree) == 1, class(spline_degree) == "integer", spline_degree >= 3L)
  stopifnot(class(return_plot_data) == "logical", length(return_plot_data) == 1, return_plot_data %in% c(TRUE, FALSE))

  # combine multiple retention times for a single sample x compound ID
  consensus_compounds <- sample_matched_compounds %>%
    dplyr::group_by(sampleId, compound_id) %>%
    dplyr::summarize(rt = stats::median(rt)) %>%
    # initialize RT deviation for a sample
    dplyr::mutate(rt_update = rt)

  rt_range <- range(sample_matched_compounds$rt)

  debugr::dwatch(msg = "Computed rt_range. [clamr<alignment.R>::estimate_rt_drift]")

  # alternate b/w: estimating a compounds true retention time (mean of RT measurements correcting for sample-to-sample drift)
  # with, fitting sample-to-sample drift from sample-level deviations from average compound RT
  for (i in 1:10) {
    debugr::dwatch(msg = paste0("Started iteration #", i, " [clamr<alignment.R>::estimate_rt_drift]"))

    # update compound rt
    compound_updates <- consensus_compounds %>%
      dplyr::group_by(compound_id) %>%
      # rt_consensus is the mean of compound retention times once the current estimate of sample-specific rt deviations has been subtracted (i.e. rt_update)
      dplyr::mutate(
        rt_consensus = stats::median(rt_update),
        rt_residual = rt - rt_consensus
      )

    debugr::dwatch(msg = paste0("compound_updates #1 computed. [clamr<alignment.R>::estimate_rt_drift]"))

    # filter compounds which are either present in very few samples or are have large RT inconsistency w.r.t rt_range

    filtered_compounds <- compound_updates %>%
      dplyr::group_by(compound_id) %>%
      dplyr::summarize(
        rt_consensus = rt_consensus[1],
        sd_rt_resid = stats::sd(rt_residual),
        sample_num = dplyr::n()
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(frac_dev = min(c(rt_consensus - rt_range[1], rt_range[2] - rt_consensus)) / diff(rt_range)) %>%
      # exclude peaks based on wildly inconsistent peaks
      dplyr::ungroup() %>%
      dplyr::filter(!(sd_rt_resid > sd_rt_resid_cutoff * diff(rt_range))) %>%
      # generate weights based on residual sd
      # threshold so precision weights don't become too massive
      dplyr::mutate(
        sd_rt_resid_thresh = pmin(pmax(sd_rt_resid, stats::quantile(.$sd_rt_resid, probs = 0.25)), stats::quantile(.$sd_rt_resid, probs = 0.75)),
        # invert weights
        compound_wt = 1 / sd_rt_resid_thresh
      ) %>%
      dplyr::select(compound_id, compound_wt)

    debugr::dwatch(msg = paste0("filtered_compounds computed. [clamr<alignment.R>::estimate_rt_drift]"))

    # for each sample, count the number of compounds that are filtered based on RT inconsistency
    sample_compound_count <- compound_updates %>%
      dplyr::semi_join(filtered_compounds, by = "compound_id") %>%
      dplyr::ungroup() %>%
      dplyr::count(sampleId) %>%
      # set up a cutoff based on # of compounds that a sample must satisfy for RT alignment
      dplyr::filter(n > spline_degree * 2)

    debugr::dwatch(msg = paste0("sample_compound_count computed. [clamr<alignment.R>::estimate_rt_drift]"))

    # some samples will not have enough points in common with any other sample to fit an update (<10 points in common with some other sample), remove these samples

    compound_updates <- compound_updates %>%
      dplyr::inner_join(filtered_compounds, by = "compound_id") %>%
      dplyr::semi_join(sample_compound_count, by = "sampleId") %>%
      dplyr::ungroup()

    debugr::dwatch(msg = paste0("compound_updates #2 computed. [clamr<alignment.R>::estimate_rt_drift]"))

    # Guard, with more detailed explanation of bug
    if (any(is.infinite(compound_updates$compound_wt))) {
      stop("[clamr<alignment.R>::estimate_rt_drift()]: some compound_updates identified with infinite weight.
           This probably means there were duplicate sample files included the set of sample files analyzed.
           Note that files may differ in name but have identical contents, which will still cause this issue.
           This is most likely an error, and that any results from this analysis should be discarded.")
    }

    # H gives a ridge penalty

    sample_gam_fits <- compound_updates %>%
      dplyr::group_by(sampleId) %>%
      dplyr::do(fit = mgcv::gam(
        formula = stats::as.formula(paste0("rt_residual ~ s(rt, k=", as.character(spline_degree), ")")),
        weights = compound_wt,
        data = .,
        H = diag(spline_ridge_penalty / spline_degree, spline_degree)
      ))

    debugr::dwatch(msg = paste0("sample_gam_fits computed. [clamr<alignment.R>::estimate_rt_drift]"))

    consensus_compounds <- update_sample_rt(
      consensus_compounds %>%
        dplyr::select(-rt_update),
      "rt", "sampleId", sample_gam_fits
    )

    debugr::dwatch(msg = paste0("consensus_compounds computed. [clamr<alignment.R>::estimate_rt_drift]"))
    debugr::dwatch(msg = paste0("Completed iteration #", i, " [clamr<alignment.R>::estimate_rt_drift]"))
  }

  debugr::dwatch(msg = "Completed 10 iterations of gam spline fit. [clamr<alignment.R>::estimate_rt_drift]")

  if (return_plot_data) {
    plot_data <- list()

    updated_consensus_compounds <- consensus_compounds %>%
      dplyr::ungroup() %>%
      dplyr::left_join(compound_updates %>%
        dplyr::select(sampleId, compound_id) %>%
        dplyr::mutate(is_retained = TRUE),
      by = c("sampleId", "compound_id")
      ) %>%
      dplyr::mutate(is_retained = ifelse(is.na(is_retained), FALSE, TRUE)) %>%
      dplyr::group_by(compound_id) %>%
      dplyr::mutate(rt_consensus = stats::median(rt_update[is_retained])) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(rt_residual = rt - rt_consensus)

    fitted_compound_deviations <- update_sample_rt(
      updated_consensus_compounds %>%
        dplyr::select(-rt_update),
      "rt", "sampleId", sample_gam_fits
    ) %>%
      dplyr::mutate(rt_deviation = rt - rt_update)

    # visualize differences between a sample RT and other features
    plot_data$rt_overall_deviations <- fitted_compound_deviations

    all_samples <- sort(unique(consensus_compounds$sampleId))
    plot_data$rt_2_sample_comparison <- expand.grid(s1 = all_samples, s2 = all_samples) %>%
      dplyr::mutate(
        s1 = ordered(s1, levels = all_samples),
        s2 = ordered(s2, levels = all_samples)
      ) %>%
      dplyr::filter(s1 < s2) %>%
      # join all features where sample is s1
      dplyr::left_join(fitted_compound_deviations %>%
        dplyr::filter(is_retained) %>%
        dplyr::mutate(sampleId = ordered(sampleId, levels = all_samples)) %>%
        dplyr::select(rt_1_raw = rt, rt_1_update = rt_update, sampleId, compound_id),
      by = c("s1" = "sampleId")
      ) %>%
      # filter to all features overlapping with sample s2
      dplyr::inner_join(fitted_compound_deviations %>%
        dplyr::filter(is_retained) %>%
        dplyr::mutate(sampleId = ordered(sampleId, levels = all_samples)) %>%
        dplyr::select(rt_2_raw = rt, rt_2_update = rt_update, sampleId, compound_id),
      by = c("s2" = "sampleId", "compound_id")
      ) %>%
      tibble::as_tibble()

    debugr::dwatch(msg = "Finished computing RT plot data. [clamr<alignment.R>::estimate_rt_drift]")
  }

  # return final RT deviation model that can be used to update RT
  # returning a key for each sample with possible starting RT value and updated RT values
  # back out retention time deviations between RT input and RT update for a dense vector of possible RT values

  output <- list()
  output$sample_gam_fits <- sample_gam_fits
  # consensus residual standard deviation between matched compounds within a full dataset
  output$rt_sd <- stats::sd(compound_updates$rt_residual[abs(compound_updates$rt_residual) < sd_rt_resid_cutoff * diff(rt_range) / 2])
  if (return_plot_data) {
    output$plot_data <- plot_data
  }

  rt_bins <- seq(0, max(sample_matched_compounds$rt) * 1.1, length.out = 1000)

  rt_bins_fit <- update_sample_rt(
    expand.grid(
      sampleId = unique(sample_matched_compounds$sampleId),
      rt = rt_bins,
      stringsAsFactors = FALSE
    ) %>%
      tibble::as_tibble(),
    "rt", "sampleId", sample_gam_fits
  )

  # cleanup by forcing RT to be non-negative and monototically increasing within sample
  rt_bins_fit <- rt_bins_fit %>%
    dplyr::arrange(sampleId, rt)

  rt_bins_fit$rt_update[rt_bins_fit$rt_update < 0] <- 0
  for (i in 2:nrow(rt_bins_fit)) {
    if (rt_bins_fit$rt_update[i] < rt_bins_fit$rt_update[i - 1] & rt_bins_fit$sampleId[i] == rt_bins_fit$sampleId[i - 1]) {
      rt_bins_fit$rt_update[i] <- rt_bins_fit$rt_update[i - 1]
    }
  }

  debugr::dwatch(msg = "Forced RT values to be non-negative and monotonically increasing within each sample. [clamr<alignment.R>::estimate_rt_drift]")

  # we thought about fitting a polynomial to this rt_update(rt) but the curve couldn't be fit exactly
  # instead, using a lookup table of rt_update(rt) on a regular grid of rt values and individual rt's can be linearly interpolated
  # from the closest two points

  # original, pre-adjustment (TODO: delete this)
  # output$rt_update_key <- rt_bins_fit

  if (rlang::is_installed("mzkitcpp")) {
    # If a previous RT alignment was performed prior to this one, map RTs back to original values.
    output$rt_update_key <- rt_bins_fit
    if (!is.null(previous_rt_update_key)) {
      updatedRts <- mzkitcpp::update_rts(previous_rt_update_key, rt_bins_fit$rt, rt_bins_fit$sampleId, debug = F)
      output$rt_update_key$rt <- updatedRts$updated_rts # the original rts (in raw files) were updated.
    }
  }

  debugr::dwatch(msg = "Updated rt_update_key using rt_bins_fit (end of function). [clamr<alignment.R>::estimate_rt_drift]")

  output
}

#' Update Sample Retention Time
#'
#' Using spline fits of observed retention time deviations from "average samples" update the retention time of a data.frame
#'
#' @param x an input data.frame (or data_frame)
#' @param rt_var the column name of an rt variable to be transformed into a new variable "rt_update"
#' @param sample_var a sample variable which matches the sample's retention time deviation gam
#' @param rt_deviation_gams a tibble of sample-specific retention time update GAMs.
update_sample_rt <- function(x, rt_var, sample_var, rt_deviation_gams) {
  stopifnot("data.frame" %in% class(x))
  stopifnot(rt_var %in% colnames(x), sample_var %in% colnames(x))
  stopifnot("tbl" %in% class(rt_deviation_gams))
  stopifnot(colnames(rt_deviation_gams) == c("sampleId", "fit"))
  if (nrow(rt_deviation_gams) != 0) {
    stopifnot(class(rt_deviation_gams$fit[[1]])[1] == "gam")
  }
  stopifnot(!(colnames(x) %in% c("rt_deviation", "rt_update")))
  if (rt_var != "rt" & "rt" %in% colnames(x)) {
    stop('rt is a reserved variable name unless it is used for "rt_var"')
  }
  if (sample_var != "sampleId" & "sampleId" %in% colnames(x)) {
    stop('sampleId is a reserved variable name unless it is used for "sample_var"')
  }

  x <- x %>%
    dplyr::rename_(
      rt = rt_var,
      sampleId = sample_var
    )

  split_x <- plyr::dlply(x, .variables = "sampleId", tibble::as_tibble)

  gam_rt_updates <- 0
  x_augment <- lapply(split_x, function(a_sample_set) {
    if (a_sample_set$sampleId[1] %in% rt_deviation_gams$sampleId) {
      # apply retention time update with GAM

      gam_rt_updates <<- gam_rt_updates + 1

      a_sample_set %>%
        dplyr::mutate(rt_deviation = mgcv::predict.gam(rt_deviation_gams$fit[rt_deviation_gams$sampleId == a_sample_set$sampleId[1]][[1]],
          newdata = a_sample_set
        )) %>%
        dplyr::mutate(rt_update = rt - rt_deviation) %>%
        dplyr::select(-rt_deviation)
    } else {
      a_sample_set %>%
        dplyr::mutate(rt_update = rt)
    }
  }) %>%
    dplyr::bind_rows()

  n_samples <- length(unique(x$sampleId))
  if (gam_rt_updates != n_samples) {
    message(round((1 - (gam_rt_updates / n_samples)) * 100), "% of samples could not be fit with a RT deviation spline - no transformation will be applied for them")
  }

  # cleanup
  x_augment %>%
    dplyr::rename_(.dots = stats::setNames(list("rt", "sampleId"), list(rt_var, sample_var)))
}

#' Estimate M/Z ppm Drift
#'
#' @description Determine the relative mass accuracy drift among experimental samples in ppm.
#'
#' @details Alternate between estimating compound-specific retention time and sample-specific retention time deviations.
#'
#' @param sample_matched_compounds tibble containing sample compound_id and rt
#' @inheritParams ms2_driven_aligner
#'
#' @return a list containing sample specific gam fits and a sample independent fit which should be used to correct sample-level overfitting
estimate_mz_ppm_drift <- function(sample_matched_compounds, return_plot_data = FALSE) {
  stopifnot(class(return_plot_data) == "logical", length(return_plot_data) == 1, return_plot_data %in% c(TRUE, FALSE))

  # combine multiple retention times for a single sample x compound ID
  updated_sample_matched_compounds <- sample_matched_compounds %>%
    dplyr::group_by(sampleId, compound_id) %>%
    dplyr::summarize(mz = stats::median(precursorMz)) %>%
    # only look at compounds where there is some disagreement about m/z
    dplyr::group_by(compound_id) %>%
    dplyr::filter(stats::sd(mz) != 0) %>%
    # initialize ppm deviation for a sample
    dplyr::mutate(ppm_deviation = 0)

  # do not perform a ppm deviation alignment on samples with < 10 shared compounds
  updated_sample_matched_compounds <- updated_sample_matched_compounds %>%
    dplyr::anti_join(updated_sample_matched_compounds %>%
      dplyr::ungroup() %>%
      dplyr::count(sampleId) %>%
      dplyr::filter(n < 10), by = "sampleId") %>%
    dplyr::ungroup()

  continue <- TRUE
  counter <- 0
  old_ppm_deviation_df <- updated_sample_matched_compounds %>%
    dplyr::ungroup() %>%
    dplyr::distinct(sampleId) %>%
    dplyr::mutate(old_ppm_deviation = 0)

  while (continue) {
    counter <- counter + 1

    # estimate consensus m/z of individual compounds
    compound_mz_consensus <- updated_sample_matched_compounds %>%
      dplyr::mutate(mz_deviation = ppm_deviation * 1e-6 * mz) %>%
      dplyr::group_by(compound_id) %>%
      dplyr::summarize(mz_compound = stats::median(mz - mz_deviation))

    # calculate ppm deviations of individual samples
    sample_ppm_deviations <- updated_sample_matched_compounds %>%
      dplyr::left_join(compound_mz_consensus, by = "compound_id") %>%
      dplyr::mutate(
        mz_deviation = mz - mz_compound,
        ppm_deviation = mz_deviation / mz * 1e6
      )

    sample_ppm_summary <- sample_ppm_deviations %>%
      dplyr::group_by(sampleId) %>%
      dplyr::summarize(ppm_deviation = stats::median(ppm_deviation))

    # check for convergence
    old_ppm_deviation_df <- sample_ppm_summary %>%
      dplyr::left_join(old_ppm_deviation_df,
        by = "sampleId"
      )

    if (sum(abs(old_ppm_deviation_df$ppm_deviation - old_ppm_deviation_df$old_ppm_deviation)) > 0.01 & counter < 20) {
      updated_sample_matched_compounds <- updated_sample_matched_compounds %>%
        dplyr::select(-ppm_deviation) %>%
        dplyr::left_join(sample_ppm_summary, by = "sampleId")

      old_ppm_deviation_df <- old_ppm_deviation_df %>%
        dplyr::select(-old_ppm_deviation) %>%
        dplyr::rename(old_ppm_deviation = ppm_deviation)
    } else {
      # convergence
      continue <- FALSE
    }
  }

  MZ_update_list <- list()
  MZ_update_list$MS1_drift <- sample_ppm_summary
  MZ_update_list$MS1_sd <- stats::sd(sample_ppm_deviations$ppm_deviation)

  if (return_plot_data) {
    MZ_update_list$plot_data <- list()
    MZ_update_list$plot_data$mz_ppm_drift <- sample_ppm_summary
    MZ_update_list$plot_data$mz_ppm_deviations <- sample_ppm_deviations
  }
  MZ_update_list
}


generate_rt_plot_data <- function(mzroll_db_con, clamr_config, ms2_groups,
                                  ms2_peak_consistency, fit_RT_update,
                                  fit_MZ_update) {
  plot_data <- list()
  # adding sample names for downstream plot labelling
  plot_data$samples <- dplyr::tbl(mzroll_db_con, "samples") %>%
    dplyr::collect()

  plot_data$sample_feature_codetection <- ms2_groups %>%
    dplyr::semi_join(ms2_peak_consistency, by = c("mz_set", "precursorCharge", "ms2_group")) %>%
    dplyr::distinct(sampleId, mz_set, ms2_group)

  plot_data <- append(plot_data, fit_RT_update$plot_data)

  if (clamr_config$MS1tol$absolute_or_relative == "relative") {
    plot_data <- append(plot_data, fit_MZ_update$plot_data)
  }

  return(plot_data)
}

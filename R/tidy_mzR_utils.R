#' Plot mzR scans
#'
#' Plot spectra for a set of scans (if more than 100 scans are provided \code{stop()} will be called).
#'
#' @param tidy_mzR output of \code{\link{tidy_mzR_from_msfile}}
#' @param scans scans to plot
#'
#' @return a ggplot2 grob
#'
#' @export
tidy_mzR_plot_scans <- function(tidy_mzR, scans) {
  tidy_mzR_test(tidy_mzR)

  stopifnot(all(class(scans) %in% c("numeric", "integer")), length(scans) > 0)
  if (length(scans) > 100) {
    stop("Too many scans provided (>100) - this function is meant to provide a detailed glimpse of a small number of scans")
  }

  reduced_scans <- tidy_mzR$scan_data %>%
    dplyr::filter(
      scan %in% scans,
      ic > 0
    )

  ggplot(reduced_scans, aes(x = mz, ymax = ic)) +
    geom_linerange(ymin = 0) +
    facet_wrap(~scan, scales = "free_y") +
    theme_minimal()
}

#' Find Major Features
#'
#' Load a mass spec file and find high intensity mzs
#'
#' @inheritParams tidy_mzR_plot_scans
#' @inheritParams test_clamr_config
#' @param ic_cutoff_fraction mzs to extract as a fraction of the maximum ic
#' @param signal_cv_cutoff mzs with a coefficient of variation of less than this value will be excluded since they are not eluting
#'
#' @return a ggplot2 grob
#'
#' @export
tidy_mzR_plot_major_features <- function(tidy_mzR, clamr_config, ic_cutoff_fraction = 1 / 20, signal_cv_cutoff = 0.2) {
  tidy_mzR_test(tidy_mzR)
  require_tolerances(clamr_config, 1L)

  stopifnot(length(ic_cutoff_fraction) == 1, class(ic_cutoff_fraction) == "numeric", ic_cutoff_fraction < 0.1, ic_cutoff_fraction > 0)
  stopifnot(length(signal_cv_cutoff) == 1, class(signal_cv_cutoff) == "numeric", signal_cv_cutoff > 0)

  # reduce to MS1 scans

  tidy_mzR$header <- tidy_mzR$header %>%
    dplyr::filter(msLevel == 1L)

  tidy_mzR$scan_data <- tidy_mzR$scan_data %>%
    dplyr::semi_join(tidy_mzR$header, by = "scan")

  # summarize distinct mzs and remove constant mzs

  fullscan_data <- tidy_mzR$scan_data %>%
    # filter for modest signal
    dplyr::filter(ic > ic_cutoff_fraction * max(.$ic) / 1000) %>%
    dplyr::arrange(mz) %>%
    # bin signals by mz tolerance
    dplyr::mutate(mz_bin = clamr::find_mz_jumps(.$mz,
      clamr_config$MS1tol$tol,
      clamr_config$MS1tol$absolute_or_relative,
      collapse = FALSE
    )) %>%
    # pool signal within a distinct mz bin
    dplyr::group_by(mz_bin, scan) %>%
    dplyr::summarize(
      mz = stats::median(mz),
      ic = sum(ic)
    ) %>%
    dplyr::group_by(mz_bin) %>%
    # generate a consensus for a distinct mz bin
    dplyr::mutate(mz = stats::median(mz)) %>%
    # remove mzs which are constant (i.e., non-eluting peaks)
    dplyr::group_by(mz_bin) %>%
    dplyr::filter(stats::sd(ic) / mean(ic) > signal_cv_cutoff) %>%
    dplyr::ungroup()

  # find major compounds

  top_mzs <- fullscan_data %>%
    dplyr::group_by(mz) %>%
    dplyr::summarize(nscans = dplyr::n(), cv = stats::sd(ic) / mean(ic), max_ic = max(ic), eic_sum = sum(ic)) %>%
    dplyr::filter(
      max_ic > ic_cutoff_fraction * max(.$max_ic),
      nscans > 20
    ) %>%
    dplyr::arrange(desc(max_ic))

  fullscan_data %>%
    dplyr::semi_join(top_mzs, by = "mz") %>%
    dplyr::left_join(tidy_mzR$header %>%
      dplyr::select(scan, retentionTime) %>%
      dplyr::mutate(retentionTime = retentionTime / 60),
    by = "scan"
    ) %>%
    dplyr::arrange(mz) %>%
    dplyr::mutate(mz = signif(mz, 6)) %>%
    dplyr::mutate(mz = factor(mz, levels = unique(.$mz))) %>%
    ggplot(aes(x = retentionTime, y = mz, height = sqrt(ic), fill = factor(as.numeric(mz) %% 3))) +
    ggjoy::geom_joy(stat = "identity", scale = 2) +
    scale_x_continuous("Retention Time", expand = c(0, 0)) +
    scale_y_discrete("M/Z") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    theme_bw()
}

#' Generated Extracted Ion Chromatographs using mzR
#'
#' Extract ion abundances per scan for a set of m/zs from a single sample.
#'
#' @inheritParams tidy_mzR_from_msfile
#' @param mz_tbl a tibble/data.frame containing one row per m/z of interest, with this mz being specified by a variable "mz"
#' @inheritParams tidy_mzR_extract_mzs
#'
#' @return mz_tbl with a nested list of EICs named "eic"
#'
#' @export
tidy_mzR_EICs <- function(ms_file_path, mzroll_db_path = NULL, mz_tbl, clamr_config) {

  # load file and extract an EIC for each standard m/z

  tidy_mzR <- tidy_mzR_from_msfile(ms_file_path,
    mzroll_db_path = mzroll_db_path
  ) # align with rt_update_key if mzroll_db_path is provided

  standard_eics <- tidy_mzR_extract_mzs(tidy_mzR,
    clamr_config,
    mz_compounds = mz_tbl$mz,
    add_missing_scans = TRUE
  )

  mz_tbl %>%
    dplyr::left_join(
      standard_eics$eic %>%
        tidyr::nest(eic = -query_mz),
      by = c("mz" = "query_mz")
    )
}

#' Extract Masses with mzR
#'
#' Open a specified file extract a compound of interest.
#'
#' @inheritParams tidy_mzR_plot_scans
#' @inheritParams test_clamr_config
#' @param mz_compounds m/z(s) to query.
#' @param add_missing_scans include scans with no signal as zeros in eic
#' @param plot_type should a plot be generated; options are: none, print and return.
#'
#' @return a list containing:
#' \itemize{
#'   \item{eic - ion counts of mz_compounds in full scan data}
#'   \item{fragmentation_events - meta information for MS2+ scans with mz_compounds as a precursorMz}
#'   \item{fragmentation_data - masses and ic of MS2+ scans}
#'   \item{summary_plot (optionally) - a summary plot will be added if \code{plot_type} = "return"}
#'   }
#'
#' @export
tidy_mzR_extract_mzs <- function(tidy_mzR, clamr_config, mz_compounds, add_missing_scans = FALSE, plot_type = "none") {
  tidy_mzR_test(tidy_mzR)
  require_tolerances(clamr_config, 1L)

  stopifnot(class(mz_compounds) == "numeric", all(mz_compounds > 0))
  stopifnot(class(add_missing_scans) == "logical", length(add_missing_scans) == 1, add_missing_scans %in% c(TRUE, FALSE))
  stopifnot(length(plot_type) == 1, class(plot_type) == "character", plot_type %in% c("none", "print", "return"))

  # EIC for mz matching this is faster than join_matching_standards()
  filter_match_fun <- switch(clamr_config$MS1tol$absolute_or_relative,
    "relative" = function(mz, mz_compound, tol) {
      abs(mz - mz_compound) < mz_compound * tol
    },
    "absolute" = function(mz, mz_compound, tol) {
      abs(mz - mz_compound) < tol
    }
  )

  ms1_scans <- tidy_mzR$scan_data %>%
    # filter for full scan events
    dplyr::filter(scan %in% tidy_mzR$header$scan[tidy_mzR$header$msLevel == 1L])

  matched_scan_data <- tibble::tibble(query_mz = mz_compounds) %>%
    dplyr::mutate(mass_matches = purrr::map(query_mz,
      tidy_mzR_extract_mzs_filter,
      ms1_scans = ms1_scans,
      clamr_config = clamr_config,
      filter_match_fun = filter_match_fun
    )) %>%
    tidyr::unnest(mass_matches)

  if (add_missing_scans) {

    # add an intensity zero entry for each scan x query_mz
    # (which will be summed with entries having signal in next step)

    matched_scan_data <- matched_scan_data %>%
      dplyr::bind_rows(
        tidy_mzR$scan_data %>%
          dplyr::filter(scan %in% tidy_mzR$header$scan[tidy_mzR$header$msLevel == 1L]) %>%
          dplyr::distinct(scan) %>%
          tidyr::crossing(query_mz = mz_compounds) %>%
          dplyr::mutate(ic = 0)
      )
  }

  eic <- matched_scan_data %>%
    # sum ic  within tolerance for each scan
    dplyr::group_by(query_mz, scan) %>%
    dplyr::summarize(ic = sum(ic)) %>%
    dplyr::left_join(tidy_mzR$header %>%
      dplyr::select(scan, rt = retentionTime) %>%
      dplyr::mutate(rt = rt / 60),
    by = "scan"
    ) %>%
    dplyr::ungroup()

  # Find all MS2+ scans matching query masses

  fragmentation_events <- tidy_mzR$header %>%
    # find MS2 events with a precursor mass within tolerance of mz_compound
    dplyr::filter(msLevel != 1L) %>%
    tidyr::crossing(query_mz = mz_compounds) %>%
    # find precursor masses within tolerance of mz_compound
    dplyr::filter(filter_match_fun(precursorMZ, query_mz, clamr_config$MS1tol$tol))

  if (plot_type != "none") {
    summary_plot <- ggplot(eic, aes(x = rt, y = ic)) +
      geom_path() +
      theme_bw() +
      geom_rug(data = fragmentation_events %>% dplyr::mutate(rt = retentionTime / 60), aes(x = rt), y = 0, sides = "b", color = "RED") +
      facet_wrap(~query_mz) +
      scale_x_continuous("Retention Time") +
      scale_y_continuous("Ion count")
  }

  output <- list(
    eic = eic,
    fragmentation_events = fragmentation_events,
    fragmentation_data = tidy_mzR$scan_data %>%
      dplyr::semi_join(fragmentation_events,
        by = "scan"
      )
  )

  if (plot_type != "none") {
    if (plot_type == "print") {
      print(summary_plot)
    } else if (plot_type == "return") {
      output$summary_plot <- summary_plot
    } else {
      stop('Invalid "plot_type"')
    }
  }

  output
}

tidy_mzR_extract_mzs_filter <- function(query_mz, ms1_scans, clamr_config, filter_match_fun) {
  ms1_scans %>%
    # find precursor masses within tolerance of mz_compound
    dplyr::filter(filter_match_fun(mz, query_mz, clamr_config$MS1tol$tol))
}

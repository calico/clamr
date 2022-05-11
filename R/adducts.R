#' Find Common Coelutions
#'
#' @description
#' Among a set of coeluting peaks, identify common mass differences as an indication of adducts or isotopologues
#'
#' @details
#' This function reads in an adduct file which contains candidate coelutions for a sample. These are generated using mzDeltas.cpp.
#'
#' @param adduct_path A file containing a set of pairs of EICs (across scan_window_length scans) with correlated intensities across scans.
#' @param coelution_n_cutoff Minimum membership size of a cluster of mass differences to retain.
#'
#' @export
find_common_coelutions <- function(adduct_path, coelution_n_cutoff = 0L) {
  stopifnot(file.exists(adduct_path))
  stopifnot(class(coelution_n_cutoff) == "integer", length(coelution_n_cutoff) == 1, coelution_n_cutoff >= 0)

  # process all distinct coelution events with a characteristic mz pair

  coelutions <- data.table::fread(
    input = adduct_path,
    col.names = c("sample", "peak_1", "peak_2", "elution_corr", "start_scan", "final_scan"),
    colClasses = c("character", "numeric", "numeric", "numeric", "integer", "integer"),
    data.table = FALSE
  ) %>%
    tibble::as_tibble() %>%
    # take the differences between peaks so that we can look for common patterns in these as indications of common adducts
    dplyr::mutate(mz_diff = peak_2 - peak_1) %>%
    dplyr::arrange(mz_diff)

  # find common mz differences between compounds
  # the grouping here is a little funky because we are matching the error of M2-M1 but M2 and M1 differ, so while error should be relative,
  # in practice we use absolute error to

  # find clearly seperable mzs

  coelutions$mz_set <- find_mz_jumps(
    sorted_mzs = coelutions$mz_diff,
    mass_accuracy_tolerance = 0.003,
    absolute_or_relative = "absolute",
    n_max_gap_inconsistencies = as.integer(ceiling(nrow(coelutions) * (1 / 1000 * 0.003))),
    collapse = FALSE
  )

  # further separate binned mass differences by density - identify local minima of KDE between peaks
  common_coelutions <- coelutions %>%
    dplyr::group_by(mz_set) %>%
    dplyr::filter(dplyr::n() >= coelution_n_cutoff) %>%
    # for each clearly seperable group find subsets corresponding to distinct peaks of mzs
    plyr::dlply(.variables = "mz_set", .fun = tibble::as_tibble) %>%
    lapply(function(a_mz_set) {
      a_mz_set$mz_subset <- find_density_minima(x = a_mz_set$mz_diff, bw = 0.001, print_plot = FALSE)
      a_mz_set
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(mz_set, mz_subset) %>%
    dplyr::filter(dplyr::n() >= coelution_n_cutoff)

  # summarize common mz differences based on mz_diff, frequency, and consistency

  coelution_mz_diffs <- common_coelutions %>%
    dplyr::group_by(mz_set, mz_subset) %>%
    dplyr::summarize(
      mean_elution_corr = mean(elution_corr),
      mz_diff_range = max(mz_diff) - min(mz_diff),
      mz_diff_median = stats::median(mz_diff),
      count = dplyr::n()
    ) %>%
    # add the mode of each group
    dplyr::inner_join(
      common_coelutions %>%
        dplyr::count(mz_set, mz_diff) %>%
        dplyr::group_by(mz_set, mz_subset) %>%
        dplyr::arrange(desc(n)) %>%
        dplyr::slice(1) %>%
        dplyr::rename(mz_diff_mode = mz_diff, mz_diff_mode_n = n),
      by = c("mz_set", "mz_subset")
    ) %>%
    dplyr::arrange(desc(count)) %>%
    dplyr::mutate(mz_diff = ifelse(mz_diff_mode_n > 10, mz_diff_mode, mz_diff_median))

  output <- list(
    coelutions = common_coelutions,
    coelution_mz_diffs = coelution_mz_diffs
  )
  class(output) <- c("list", "coelution")
  output
}

#' Identify Coelutions
#'
#' Label common mass differences using either a library of commonly seen adducts and isotopologues or optimize chemical composition as a mixed integer linear programming problem.
#'
#' @param coelutions A list containing coelutions and coelution_mz_diffs of the \code{coelution} class produced by \code{find_common_coelutions}
#' @param type
#' \itemize{
#'  \item{\code{library}: match mzs to a list of knowns based on a specified mass accuracy tolerance -- required objects coelution_library, amu_mass_accuracy_tolerance}
#'  \item{\code{denovo}: match each mz to linear combinations of isotopes -- required objects: isotope_summary}
#' }
#' @param coelution_n_cutoff Minimum membership size of a mz difference group to be tested for IDs
#' @param quietly suppress print statements
#' @param coelution_library a data.frame containing the following required variables
#' \itemize{
#'  \item{\code{mz}: numeric mz differences in amu}
#'  \item{\code{name}: character}
#'  \item{\code{type}: the type of coelution event -- e.g., isotopologue, adduct, distinct_compound}
#'  }
#' @param amu_tolerance search tolerance for matches in amu -- will report the highest priority match and otherwise the closest match
#' @inheritParams estimate_chemical_formula
#' @param ... additional arguments to pass to coelution identification
#'
#' @export
identify_coelutions <- function(coelutions, type = "denovo", coelution_n_cutoff = 0L, quietly = FALSE, coelution_library = NULL,
                                amu_tolerance = NULL, isotope_summary = NULL, ...) {

  # test inputs
  if (!("coelution" %in% class(coelutions))) {
    stop('"coelutions" must contain the "coelution" class')
  }

  stopifnot(class(type) == "character", length(type) == 1)
  type_options <- c("library", "denovo")
  if (!(type %in% type_options)) {
    stop(paste0(type, ' is not a valid option for "type"; valid options are: ', paste(type_options, collapse = ", ")))
  } else if (type == "library") {
    stopifnot(class(coelution_library) != "NULL", "data.frame" %in% class(coelution_library), all(c("mz_diff", "name", "type") %in% colnames(coelution_library)))
    stopifnot(class(amu_tolerance) == "numeric", length(amu_tolerance) == 1, amu_tolerance > 0)
  } else if (type == "denovo") {
    stopifnot(class(isotope_summary) != "NULL", "data.frame" %in% class(isotope_summary), all(c("label", "amu") %in% colnames(isotope_summary)))
  } else {
    stop(paste0(type, " is in type_options but not contain any type-specific behaviors"))
  }

  stopifnot(class(coelution_n_cutoff) == "integer", length(coelution_n_cutoff) == 1, coelution_n_cutoff >= 0)
  stopifnot(class(quietly) == "logical", length(quietly) == 1, quietly %in% c(TRUE, FALSE))

  # capture optional inputs
  dots <- list(...)

  # filter to a set of query masses based on coelution_n_cutoff

  query_coelutions <- coelutions$coelution_mz_diffs %>%
    dplyr::filter(count >= coelution_n_cutoff)

  if (nrow(query_coelutions) == 0) {
    warning(paste0('no enries were present more than "coelution_n_cutoff" times (', coelution_n_cutoff, ")"))
    return(
      # return an empty tibble with appropriate classes
      tibble::tibble(
        mz_set = 1L, mz_subset = 1L, mean_elution_corr = 0.5,
        mz_diff_range = 0.5, mz_diff_median = 0.5, count = 1L,
        mz_diff_mode = 0.5, mz_diff_mode_n = 1L,
        query_mz_diff = 0.5, name = "", formula = "", type = "",
        mode_specific = "", analyte_specific = TRUE,
        standard_mz_diff = 0.5, is_loss = TRUE,
        match_distance = 0.5
      ) %>%
        dplyr::slice(-1)
    )
  } else if (!quietly) {
    message(paste0(nrow(query_coelutions), " mass differences will be searched"))
  }

  if (type == "library") {
    # search a library of annotated mass differences
    join_matching_standards(
      query_dataset = query_coelutions,
      standard_dataset = coelution_library,
      variable_tolerances = tibble::tibble(variable = "mz_diff", tolerance = amu_tolerance, relative_or_absolute = "absolute"),
      threshold = 1
    )
  } else if (type == "denovo") {
    # search for a molecular formula based on linear combinations of allowable isotopes

    search_options <- dots[intersect(names(formals(estimate_chemical_formula)), names(dots))]

    adduct_formulas <- parallel::mclapply(query_coelutions$mz_diff, function(a_mass) {
      do.call(estimate_chemical_formula, append(list(query_mass = a_mass, isotope_summary = isotope_summary), search_options))
    }, mc.cores = pmax(1, parallel::detectCores() - 2)) %>%
      purrr::transpose() %>%
      purrr::map(dplyr::bind_rows)

    adduct_formulas$objective
  }
}

#' Apply Coelution Labels
#'
#' @inheritParams identify_coelutions
#' @param coelution_labels a data_frame containing:
#' \itemize{
#'   \item{mz_set}
#'   \item{mz_subset}
#'   \item{library_formula - a formula label}
#'   \item{library_is_loss - TRUE/FALSE - whether the coelution difference is due to a gain or loss of the monoisotopic / non-adduct form}
#'   \item{library_type - label to apply to the connections - valid labels are "adduct", "isotopologue" or "distinct compound"}
#'   }
#'
#' @return a coelution object with labelled mz differences (added variables: library_formula, library_is_loss, library_type)
#'
#' @export
apply_coelution_labels <- function(coelutions, coelution_labels) {

  # test inputs
  if (!("coelution" %in% class(coelutions))) {
    stop('"coelutions" must contain the "coelution" class')
  }

  stopifnot(all(colnames(coelution_labels) == c("mz_set", "mz_subset", "library_formula", "library_is_loss", "library_type")))

  nonsupported_types <- setdiff(unique(coelution_labels$library_type), c("adduct", "isotopologue", "distinct compound"))
  if (length(nonsupported_types) != 0) {
    warning(paste0(paste(nonsupported_types, collapse = ", "), " are not supported in library_type"))
  }

  coelutions$coelution_mz_diffs <- coelutions$coelution_mz_diffs %>%
    dplyr::inner_join(coelution_labels, by = c("mz_set", "mz_subset")) %>%
    dplyr::ungroup()

  reduced_coelutions <- coelutions$coelutions %>%
    dplyr::inner_join(coelutions$coelution_mz_diffs %>%
      dplyr::select(mz_set, mz_subset, library_formula, library_is_loss, library_type),
    by = c("mz_set", "mz_subset")
    ) %>%
    dplyr::ungroup()

  stopifnot(all(reduced_coelutions$library_is_loss %in% c(TRUE, FALSE)))

  # flip peak_1 and peak_2 is library_is_loss is TRUE, indicating that the monoisotopic / non-adduct mass is the larger mass (normally it is smaller)
  coelutions$coelutions <- reduced_coelutions %>%
    dplyr::filter(!library_is_loss) %>%
    dplyr::bind_rows(
      reduced_coelutions %>%
        dplyr::filter(library_is_loss) %>%
        dplyr::rename(peak_2 = peak_1, peak_1 = peak_2)
    )

  coelutions
}

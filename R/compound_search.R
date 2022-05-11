#' Nest Peak Groups
#'
#' Combine the peak and scan data for each peakgroup into a 1-row per peakgroup table containing nested peak and scan lists.
#'
#' @inheritParams test_mzroll_db_con_schema
#' @inheritParams test_clamr_config
#'
#' @return a tibble containing nested peaks and scans for each peak group
#'
#' @export
nest_peakgroup_features <- function(mzroll_db_con) {
  test_mzroll_db_con_schema(mzroll_db_con)

  has_parentGroupId_tbl <- dplyr::tbl(
    mzroll_db_con,
    dbplyr::sql("SELECT COUNT(*) AS hasParentGroupIdColumn FROM pragma_table_info('peakgroups') WHERE name='parentGroupId'")
  ) %>%
    dplyr::collect()

  if (has_parentGroupId_tbl$hasParentGroupIdColumn == 0) {
    stop("this mzrollDB does not contain the column 'peakGroupId' in the 'peakgroups' table, indicating that it has already been mutated, most likely by a previous mzkit pipeline_standard_search analysis.  A clamr_db file that has already been searched via pipeline_standard_search cannot be searched again.")
  }

  peakgroups_data <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT groupId, parentGroupId, adductName FROM peakgroups")) %>%
    dplyr::mutate(adductName = ifelse(adductName == "[M-H]", "[M-H]-", adductName)) %>%
    dplyr::mutate(adductName = ifelse(adductName == "[M+H]", "[M+H]+", adductName)) %>%
    dplyr::collect()

  # match clamr peaks
  peak_data <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT peakId, groupId, sampleId, rt as peakRt, rtmin, rtmax, peakMz, mzmin, mzmax, minscan, maxscan, peakAreaTop, quality FROM peaks")) %>%
    dplyr::collect() %>%
    dplyr::semi_join(peakgroups_data, by = "groupId")

  scan_data <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT sampleId, scan, rt, precursorMz, precursorCharge, precursorPurity, data FROM scans WHERE mslevel == 2")) %>%
    dplyr::collect()

  # join peaks to scans when there is a scan range and mz range match for a given sample
  peak_scan_matches <- tidyr::nest_legacy(peak_data %>%
    dplyr::select(peakId, groupId, sampleId, peakMz, mzmin, mzmax, minscan, maxscan), -sampleId, .key = "peaks") %>%
    dplyr::left_join(tidyr::nest_legacy(scan_data, -sampleId, .key = "scans"), by = "sampleId") %>%
    dplyr::mutate(matches = purrr::map2(peaks, scans, mz_scan_join_fxn)) %>%
    tidyr::unnest_legacy(matches) %>%
    # require that each scan was matched to a peak
    dplyr::filter(!is.na(peakId))

  nested_peakgroup_features <- peakgroups_data %>%
    dplyr::left_join(tidyr::nest_legacy(peak_data, -groupId, .key = "peaks"), by = "groupId") %>%
    dplyr::left_join(tidyr::nest_legacy(peak_scan_matches, -groupId, .key = "scans"), by = "groupId")

  nested_peakgroup_features
}

#' Join Matching Standards
#'
#' @description
#' Find all entries in a standard dataset which match features in a query dastaset based on the consistency of 1 or more variables.
#'
#' @details
#' Distance is calculated on variables that have been rescaled to their tolerances (variable/tolerance) so that deviations in each dimension are equivalent.
#'
#' @param query_dataset containing candidate features with 1+ match variables.
#' @param standard_dataset containing standards with 1+ variables which match \code{query_dataset}.
#' @param variable_tolerances a three column tibble containing "variable", "tolerance" and "relative_or_absolute"
#' @param distance_method distance measure between scaled variables (x/tol): e.g., manhattan, euclidean
#' @param threshold distance cutoff to report a match
#'
#' @return an inner_join tibble of the query_dataset and standard_dataset
#'
#' @examples
#' library(dplyr)
#'
#' query_dataset <- tibble::tibble(query_id = 1:3, amu = 1:3)
#' standard_dataset <- clamr::isotope_summaries %>% dplyr::select(z, label, amu)
#' variable_tolerances <- tibble::tibble(variable = "amu", tolerance = 0.01, relative_or_absolute = "absolute")
#' join_matching_standards(query_dataset, standard_dataset, variable_tolerances, distance_method = "manhattan", threshold = 1)
#'
#' query_dataset <- tibble::tibble(query_id = 1:3, amu = 1:3, rt = c(10, 15, 3))
#' standard_dataset <- clamr::isotope_summaries %>%
#'   dplyr::select(z, label, amu) %>%
#'   dplyr::mutate(rt = 1:n())
#' variable_tolerances <- tibble::tribble(
#'   ~variable, ~tolerance, ~relative_or_absolute,
#'   "amu", 0.01, "absolute",
#'   "rt", 5, "absolute"
#' )
#' join_matching_standards(query_dataset, standard_dataset, variable_tolerances, distance_method = "manhattan", threshold = 4)
#'
#' # matching with ppm
#' query_dataset <- tibble::tibble(query_id = 1:4, mz = c(1, 10, 100, 1000))
#' standard_dataset <- tibble::tibble(library_id = 1:4, mz = c(1.001, 10.001, 100.001, 1000.001))
#' variable_tolerances <- tibble::tibble(variable = "mz", tolerance = 10e-6, relative_or_absolute = "relative")
#' join_matching_standards(query_dataset, standard_dataset, variable_tolerances, distance_method = "manhattan", threshold = 10)
#'
#' query_dataset <- tibble::tibble(query_id = 1:4, mz = c(1, 10, 100, 1000), rt = 1:4)
#' standard_dataset <- tibble::tibble(library_id = 1:4, mz = c(1.001, 10.001, 100.001, 1000.001), rt = 1:4 + 0.5)
#' variable_tolerances <- tibble::tribble(
#'   ~variable, ~tolerance, ~relative_or_absolute,
#'   "mz", 10e-6, "relative",
#'   "rt", 1, "absolute"
#' )
#' join_matching_standards(query_dataset, standard_dataset, variable_tolerances, distance_method = "manhattan", threshold = 10)
#'
#' @export
join_matching_standards <- function(query_dataset, standard_dataset, variable_tolerances, distance_method = "manhattan", threshold = nrow(variable_tolerances)) {
  checkmate::assertDataFrame(query_dataset)
  checkmate::assertDataFrame(standard_dataset)
  checkmate::assertDataFrame(variable_tolerances)
  stopifnot(variable_tolerances$relative_or_absolute %in% c("relative", "absolute"))

  if (!all(variable_tolerances$variable %in% colnames(query_dataset))) {
    stop(paste(setdiff(variable_tolerances$variable, colnames(query_dataset)), collapse = ", "), ' were not found in query_dataset. "variable_tolerances" should be defined for all match variables.')
  }
  if (!all(variable_tolerances$variable %in% colnames(query_dataset))) {
    stop(paste(setdiff(variable_tolerances$variable, colnames(standard_dataset)), collapse = ", "), ' were not found in standard_dataset. "variable_tolerances" should be defined for all match variables.')
  }

  query_dataset <- query_dataset %>%
    dplyr::ungroup()
  minimal_queries <- query_dataset %>%
    dplyr::ungroup() %>%
    dplyr::mutate(.query_row = 1:dplyr::n()) %>%
    dplyr::select(!!!rlang::syms(c(".query_row", variable_tolerances$variable))) %>%
    tidyr::drop_na()

  standard_dataset <- standard_dataset %>%
    dplyr::ungroup()
  minimal_library <- standard_dataset %>%
    dplyr::mutate(.standard_row = 1:dplyr::n()) %>%
    dplyr::select(!!!rlang::syms(c(".standard_row", variable_tolerances$variable))) %>%
    tidyr::drop_na()

  # find matches between

  if (nrow(variable_tolerances) == 1) {
    # if there is only 1 variable then just calculate the scaled difference between query and library

    vectorized_scaled_diff <- function(x, y, tolerance, relative_or_absolute = "absolute") {
      if (relative_or_absolute == "absolute") {
        abs(x - y) / tolerance
      } else if (relative_or_absolute == "relative") {
        abs(x - y) / (tolerance * x)
      } else {
        stop("invalid value of relative_or_absolute")
      }
    }

    # fuzzyjoin gives a warning about ... which doesn't seem to affect anything
    # look like how fuzzyjoin is calling an underlying function

    library_matches <- fuzzyjoin::fuzzy_join(
      x = minimal_queries,
      y = minimal_library,
      by = variable_tolerances$variable,
      match_fun = function(x, y, tolerance, relative_or_absolute, threshold) {
        vectorized_scaled_diff(x, y, tolerance, relative_or_absolute) < threshold
      },
      threshold = threshold,
      tolerance = variable_tolerances$tolerance,
      relative_or_absolute = variable_tolerances$relative_or_absolute
    )
    if (nrow(library_matches) == 0) {
      return(NULL)
    } else {
      mutate_call <- lazyeval::interp(~ vectorized_scaled_diff(x, y, tolerance, relative_or_absolute),
        x = as.name(paste0(variable_tolerances$variable, ".x")),
        y = as.name(paste0(variable_tolerances$variable, ".y")),
        tolerance = variable_tolerances$tolerance,
        relative_or_absolute = variable_tolerances$relative_or_absolute
      )

      library_matches <- library_matches %>%
        dplyr::mutate_(.dots = setNames(list(mutate_call), "match_distance"))
    }
  } else {
    # join by multiple columns
    distance_function <- switch(distance_method,
      manhattan = function(D) {
        rowSums(abs(D))
      },
      euclidean = function(D) {
        sqrt(rowSums(D^2))
      },
      maximum = function(D) {
        apply(abs(D), 1, max)
      },
      stop(paste0(distance_method, ' is not a defined "distance_method" when only one match variable exists'))
    )

    scaled_distance_match <- function(X, Y) {
      x_y_diffs <- lapply(1:nrow(variable_tolerances), function(j) {
        x <- as.numeric(X[, j] %>% unlist() %>% unname())
        y <- as.numeric(Y[, j] %>% unlist() %>% unname())

        if (variable_tolerances$relative_or_absolute[variable_tolerances$variable == colnames(X)[j]] == "absolute") {
          (x - y) / variable_tolerances$tolerance[variable_tolerances$variable == colnames(X)[j]]
        } else {
          (x - y) / (variable_tolerances$tolerance[variable_tolerances$variable == colnames(X)[j]] * x)
        }
      }) %>%
        unlist() %>%
        matrix(ncol = nrow(variable_tolerances), byrow = FALSE) %>%
        distance_function()

      tibble::tibble(x_y_diffs < threshold, match_distance = x_y_diffs)
    }

    library_matches <- fuzzyjoin::fuzzy_join(minimal_queries,
      minimal_library,
      multi_by = variable_tolerances$variable,
      multi_match_fun = scaled_distance_match,
      distance_function = distance_function,
      threshold = threshold
    )
  }

  if (nrow(library_matches) != 0) {
    query_dataset %>%
      dplyr::rename_(.dots = as.list(setNames(variable_tolerances$variable, paste0("query_", variable_tolerances$variable)))) %>%
      dplyr::slice(library_matches$.query_row) %>%
      dplyr::bind_cols(
        standard_dataset %>%
          dplyr::rename_(.dots = as.list(setNames(variable_tolerances$variable, paste0("standard_", variable_tolerances$variable)))) %>%
          dplyr::slice(library_matches$.standard_row)
      ) %>%
      dplyr::mutate(match_distance = library_matches$match_distance)
  } else {
    NULL
  }
}

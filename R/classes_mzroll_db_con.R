#' mzrollDB file (SQLite)
#'
#' Connect to an mzrollDB SQLite database.
#'
#' @param mzroll_db_path path to a .mzrollDB SQLite file.
#' @inheritParams create_sqlite_con
#'
#' @return \code{mzroll_db_con}: a connection to the mzroll/sqlite database
#'
#' @export
mzroll_db_sqlite <- function(mzroll_db_path, unlock = TRUE) {
  stopifnot(class(mzroll_db_path) == "character", length(mzroll_db_path) == 1)
  stopifnot(class(unlock) == "logical", length(unlock) == 1, unlock %in% c(TRUE, FALSE))

  if (!file.exists(mzroll_db_path)) {
    stop("mzroll_db_path: ", mzroll_db_path, " not found")
  }

  mzroll_db_con <- create_sqlite_con(mzroll_db_path, unlock = unlock)
  test_mzroll_db_con_schema(mzroll_db_con)

  return(mzroll_db_con)
}

#' Test mzrollDB Connection Schema
#'
#' @param mzroll_db_con a connection to a mzroll database as produced by \code{\link{mzroll_db_sqlite}}
#'
#' @export
test_mzroll_db_con_schema <- function(mzroll_db_con) {
  if (!file.exists(mzroll_db_con@dbname)) {
    stop(paste0("database path: ", mzroll_db_con@dbname, " does not exist"))
  }

  # test for required tables
  db_tables <- DBI::dbListTables(mzroll_db_con)
  standard_tables <- c("peaks", "samples", "scans")
  if (!all(standard_tables %in% db_tables)) {
    stop("required tables: ", paste(setdiff(standard_tables, db_tables), collapse = ", "), " were not found in mzroll_db_con")
  }

  # summarize the schema
  mz_roll_db_schema <- sql_get_schema(mzroll_db_con)

  # test the schema for appropriate classes
  validate_mzroll_db_schema(mz_roll_db_schema$db_head)

  # test more specific formatting
  if (!all(grepl("^(\\[[0-9.]+,[0-9.e+-]+\\])+$", mz_roll_db_schema$db_head$scans$data))) {
    stop("scan$data is an unexpected format; an example of the required format is: [64.9774475,63309][256.905396,22314][90.9764404,12070]")
  }

  # test for uniqueness of groupIds and sampleIds
  test_mzroll_db_con_schema_pks(mzroll_db_con)

  invisible(0)
}

#' Validate mzrollDB Schema
#'
#' Compare a list of tables in mzroll_data_list to the mzrollDB schema
#' to determine whether the list is compatible with the list schema
#'
#' @param mzroll_data_list a list of tables matching the mzrollDB schema or a subset of the required tables.
#' @param tables either "all" to compare each table to the schema and ensure that all expected tables are present or
#'   a character vector of tables to compare.
#'
#' @return 0 invisibly (used for side-effects)
validate_mzroll_db_schema <- function(mzroll_data_list, tables = "all") {
  stopifnot("list" %in% class(mzroll_data_list))
  stopifnot(class(tables) == "character")

  mz_roll_var_classes <- lapply(seq_along(mzroll_data_list), function(table_i) {
    mzroll_data_list[[table_i]] %>%
      purrr::map(class) %>%
      {
        tibble::tibble(table = names(mzroll_data_list)[table_i], variable = names(.), class = sapply(., function(x) {
          paste(x, collapse = "-")
        }))
      }
  }) %>%
    dplyr::bind_rows()

  # test for required variables and matching classes

  mz_roll_schema_classes <- mzroll_schema_required_classes()

  if (length(tables) != 1 || tables != "all") {

    # define tables to evaluate
    undefined_tables <- setdiff(tables, mz_roll_schema_classes$table)
    if (length(undefined_tables) != 0) {
      stop("tables to be tested are not specified in mzrollDB schema: ", paste(undefined_tables, collapse = ", "))
    }

    mz_roll_schema_classes <- mz_roll_schema_classes %>%
      dplyr::filter(table %in% tables)
  }

  required_variable_table <- mz_roll_schema_classes %>%
    dplyr::filter(variable_is_required) %>%
    dplyr::left_join(mz_roll_var_classes, by = c("table", "variable"))

  required_vars_missing <- required_variable_table %>%
    dplyr::filter(is.na(class)) %>%
    dplyr::mutate(violation = "required variable is missing")

  vars_inconsistent_classes <- mz_roll_var_classes %>%
    dplyr::inner_join(mz_roll_schema_classes, by = c("table", "variable")) %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::filter(class != required_class) %>%
    dplyr::mutate(violation = dplyr::case_when(
      variable_is_required == TRUE ~ "required variable - variable class does not match required class",
      variable_is_required == FALSE ~ "optional variable - variable class does not match required class"
    ))

  if (nrow(required_vars_missing) != 0 | nrow(vars_inconsistent_classes) != 0) {
    schema_problems <- required_vars_missing %>%
      dplyr::bind_rows(vars_inconsistent_classes)

    schema_warnings <- schema_problems %>%
      dplyr::filter(!variable_is_required)

    if (nrow(schema_warnings) != 0) {
      print(schema_warnings)
      warning(nrow(schema_warnings), " schema warnings\n")
    }

    schema_errors <- schema_problems %>%
      dplyr::filter(variable_is_required)

    if (nrow(schema_errors) != 0) {
      print(schema_errors)
      stop(nrow(schema_errors), " schema errors\n")
    }
  }

  invisible(0)
}


#' CLaM R Summary
#'
#' Summarize the state of an mzrollDB dataset
#'
#' @param mzroll_db_con Test for the appropriateness of mzroll_db_con's sql schema
#' @param clamr_config a list containing parameters mass tolerances and other dataset parameters
#'
#' @export
clamr_summary <- function(mzroll_db_con, clamr_config) {
  db_tables <- DBI::dbListTables(mzroll_db_con)

  ms_file_summary <- tibble::frame_data(~Metric, ~Value)

  if ("peakgroups" %in% db_tables) {
    ms_file_summary <- ms_file_summary %>%
      dplyr::bind_rows(tibble::frame_data(
        ~Metric, ~Value,
        "Peak Groups", "",
        "number of peak groups", dplyr::tbl(mzroll_db_con, dplyr::sql("SELECT Count(*) FROM peakgroups")) %>% dplyr::collect() %>% unlist() %>% unname() %>%
          {
            format(., big.mark = ",", scientific = FALSE)
          },
        "number of compound IDs", dplyr::tbl(mzroll_db_con, dplyr::sql("SELECT Count(*) FROM peakgroups WHERE compoundName != ''")) %>% dplyr::collect() %>% unlist() %>% unname() %>%
          {
            format(., big.mark = ",", scientific = FALSE)
          },
        "number of adduct IDs", dplyr::tbl(mzroll_db_con, dplyr::sql("SELECT Count(*) FROM peakgroups WHERE adductName != ''")) %>% dplyr::collect() %>% unlist() %>% unname() %>%
          {
            format(., big.mark = ",", scientific = FALSE)
          },
        "", ""
      ))
  }

  if ("peaks" %in% db_tables) {
    ms_file_summary <- ms_file_summary %>%
      dplyr::bind_rows(tibble::frame_data(
        ~Metric, ~Value,
        "Peaks", "",
        "MZ range", dplyr::tbl(mzroll_db_con, dplyr::sql("SELECT MIN(mzmin), MAX(mzmax) FROM peaks")) %>% dplyr::collect() %>% unlist() %>% unname() %>% paste(., collapse = "-"),
        "RT range", dplyr::tbl(mzroll_db_con, dplyr::sql("SELECT MIN(rtmin), MAX(rtmax) FROM peaks")) %>% dplyr::collect() %>% unlist() %>% unname() %>% paste(., collapse = "-"),
        "", ""
      ))
  }

  if ("samples" %in% db_tables) {
    all_samples <- dplyr::tbl(mzroll_db_con, dplyr::sql("SELECT name FROM samples")) %>%
      dplyr::collect() %>%
      unlist() %>%
      unname()

    ms_file_summary <- ms_file_summary %>%
      dplyr::bind_rows(tibble::frame_data(
        ~Metric, ~Value,
        "Samples", "",
        "number of samples", length(all_samples),
        "number of blanks", sum(stringr::str_detect(tolower(all_samples), "blank")),
        "number of controls", sum(stringr::str_detect(tolower(all_samples), "ctl|control")),
        "number of standards", sum(stringr::str_detect(tolower(all_samples), "std|standard")),
        "", ""
      ))
  }

  if ("scans" %in% db_tables) {
    MSN_counts <- dplyr::tbl(mzroll_db_con, dplyr::sql("SELECT mslevel, Count(*) FROM scans GROUP BY mslevel")) %>% dplyr::collect()
    colnames(MSN_counts) <- c("Metric", "Value")
    MSN_counts$Metric <- paste("MS", MSN_counts$Metric, " scans", sep = "")
    MSN_counts$Value <- format(MSN_counts$Value, big.mark = ",", scientific = FALSE)

    ms_file_summary <- ms_file_summary %>%
      dplyr::bind_rows(tibble::frame_data(
        ~Metric, ~Value,
        "Scans", ""
      )) %>%
      dplyr::bind_rows(MSN_counts) %>%
      dplyr::bind_rows(tibble::frame_data(
        ~Metric, ~Value,
        "", ""
      ))
  }

  acquisition_par_summary <- lapply(seq_along(clamr_config), function(i) {
    a_parameter <- clamr_config[[i]]
    a_parameter_name <- names(clamr_config)[i]

    if (length(a_parameter) == 1) {
      tibble::tibble(Metric = a_parameter_name, Value = a_parameter)
    } else {
      tibble::tibble(
        Metric = a_parameter_name,
        Value = unlist(a_parameter) %>%
          {
            paste(names(.), ": ", unname(.), sep = "")
          } %>%
          {
            paste(., collapse = ", ")
          }
      )
    }
  }) %>%
    dplyr::bind_rows()

  output <- ms_file_summary %>%
    dplyr::bind_rows(tibble::tibble(Metric = "Acquisition Params", Value = "")) %>%
    dplyr::bind_rows(acquisition_par_summary)

  print(output, n = 100)
}

mzroll_schema_required_classes <- function() {
  tibble::frame_data(
    ~table, ~variable, ~required_class, ~variable_is_required,
    "peakgroups", "groupId", "integer", T,
    "peakgroups", "parentGroupId", "integer", T,
    "peakgroups", "tagString", "character", T,
    "peakgroups", "metaGroupId", "integer", T,
    "peakgroups", "expectedRtDiff", "numeric", F,
    "peakgroups", "groupRank", "numeric", F,
    "peakgroups", "label", "character", F,
    "peakgroups", "type", "integer", F,
    "peakgroups", "srmId", "character", F,
    "peakgroups", "ms2EventCount", "integer", F,
    "peakgroups", "ms2Score", "numeric", F,
    "peakgroups", "adductName", "character", F,
    "peakgroups", "compoundId", "character", F,
    "peakgroups", "compoundName", "character", F,
    "peakgroups", "compoundDB", "character", F,
    "peakgroups", "searchTableName", "character", F,
    "peaks", "peakId", "integer", T,
    "peaks", "groupId", "integer", T,
    "peaks", "sampleId", "integer", T,
    "peaks", "pos", "integer", T,
    "peaks", "minpos", "integer", T,
    "peaks", "maxpos", "integer", T,
    "peaks", "rt", "numeric", T,
    "peaks", "rtmin", "numeric", T,
    "peaks", "rtmax", "numeric", T,
    "peaks", "mzmin", "numeric", T,
    "peaks", "mzmax", "numeric", T,
    "peaks", "scan", "integer", T,
    "peaks", "minscan", "integer", T,
    "peaks", "maxscan", "integer", T,
    "peaks", "peakArea", "numeric", F,
    "peaks", "peakAreaCorrected", "numeric", F,
    "peaks", "peakAreaTop", "numeric", T,
    "peaks", "peakAreaFractional", "numeric", F,
    "peaks", "peakRank", "numeric", F,
    "peaks", "peakIntensity", "numeric", F,
    "peaks", "peakBaseLineLevel", "numeric", F,
    "peaks", "peakMz", "numeric", T,
    "peaks", "medianMz", "numeric", T,
    "peaks", "baseMz", "numeric", T,
    "peaks", "quality", "numeric", T,
    "peaks", "width", "integer", F,
    "peaks", "gaussFitSigma", "numeric", F,
    "peaks", "gaussFitR2", "numeric", F,
    "peaks", "noNoiseObs", "integer", F,
    "peaks", "noNoiseFraction", "numeric", F,
    "peaks", "symmetry", "numeric", F,
    "peaks", "signalBaselineRatio", "numeric", F,
    "peaks", "groupOverlap", "numeric", F,
    "peaks", "groupOverlapFrac", "numeric", F,
    "peaks", "localMaxFlag", "numeric", F,
    "peaks", "fromBlankSample", "integer", F,
    "peaks", "label", "integer", F,
    "samples", "sampleId", "integer", T,
    "samples", "name", "character", T,
    "samples", "filename", "character", T,
    "samples", "setName", "character", T,
    "scans", "id", "integer", T,
    "scans", "sampleId", "integer", T,
    "scans", "scan", "integer", T,
    "scans", "fileSeekStart", "integer", F,
    "scans", "fileSeekEnd", "integer", F,
    "scans", "mslevel", "integer", T,
    "scans", "rt", "numeric", T,
    "scans", "precursorMz", "numeric", T,
    "scans", "precursorCharge", "integer", T,
    "scans", "precursorIc", "numeric", T,
    "scans", "minmz", "numeric", T,
    "scans", "maxmz", "numeric", T,
    "scans", "data", "character", T
  )
}

#' Expand With Mzroll Defaults
#'
#' @param mzroll_data a list of tables corresponding to an mzrollDB SQL
#'   schema, but possibly missing some variables.
#'
#' @return mzroll_data with default values added for variables which were
#'   not initially included
#'
#' @examples
#' library(dplyr)
#'
#' peakgroups <- tibble::tibble(groupId = 1:10, compoundName = LETTERS[1:10])
#' samples <- tibble::tibble(sampleId = 1:10)
#' peaks <- tidyr::crossing(peakgroups, samples) %>%
#'   dplyr::mutate(peakId = 1:dplyr::n())
#'
#' mzroll_data <- list(
#'   peakgroups = peakgroups,
#'   samples = samples,
#'   peaks = peaks
#' )
#'
#' expand_with_mzroll_defaults(mzroll_data)
#'
#' @export
expand_with_mzroll_defaults <- function(mzroll_data) {
  checkmate::assertClass(mzroll_data, "list")
  stopifnot(length(names(mzroll_data)) == length(mzroll_data))

  # verify that all tables have defaults
  mzroll_defaults <- mzroll_defaults()
  purrr::walk(names(mzroll_data), function(x) {
    checkmate::assertClass(mzroll_data[[x]], "data.frame")
    checkmate::assertChoice(x, names(mzroll_defaults))
  })

  augmented_mzroll_data <- purrr::map(names(mzroll_data), function(x) {
    var_order <- dplyr::filter(
      mzroll_schema_required_classes(),
      table == x
    )[["variable"]]

    missing_null_fields <- dplyr::select(
      mzroll_defaults[[x]],
      !!!rlang::syms(setdiff(
        colnames(mzroll_defaults[[x]]),
        colnames(mzroll_data[[x]])
      ))
    )

    tidyr::crossing(
      mzroll_data[[x]],
      missing_null_fields
    ) %>%
      dplyr::select(!!!rlang::syms(var_order))
  })

  names(augmented_mzroll_data) <- names(mzroll_data)

  # check that all required fields are present
  validate_mzroll_db_schema(
    augmented_mzroll_data,
    tables = names(augmented_mzroll_data)
  )

  return(augmented_mzroll_data)
}

#' MzRoll Defaults
#'
#' Default values for mzroll tables when defaults are possible. groupId,
#'   sampleId, and peakId do not have defaults.
mzroll_defaults <- function() {
  peakgroups <- tibble::tibble(
    parentGroupId = 0L,
    tagString = "",
    metaGroupId = 0L,
    expectedRtDiff = -1,
    groupRank = 0,
    label = "",
    type = NA_integer_,
    srmId = "",
    ms2EventCount = NA_integer_,
    ms2Score = 0,
    adductName = NA_character_,
    compoundId = NA_character_,
    compoundName = NA_character_,
    compoundDB = NA_character_,
    searchTableName = NA_character_
  )

  samples <- tibble::tibble(
    name = NA_character_,
    filename = NA_character_,
    setName = NA_character_
  )

  peaks <- tibble::tibble(
    pos = NA_integer_,
    minpos = NA_integer_,
    maxpos = NA_integer_,
    rt = NA_real_,
    rtmin = NA_real_,
    rtmax = NA_real_,
    mzmin = NA_real_,
    mzmax = NA_real_,
    scan = NA_integer_,
    minscan = NA_integer_,
    maxscan = NA_integer_,
    peakArea = NA_real_,
    peakAreaCorrected = NA_real_,
    peakAreaTop = NA_real_,
    peakAreaFractional = NA_real_,
    peakRank = NA_real_,
    peakIntensity = NA_real_,
    peakBaseLineLevel = NA_real_,
    peakMz = NA_real_,
    medianMz = NA_real_,
    baseMz = NA_real_,
    quality = NA_real_,
    width = NA_integer_,
    gaussFitSigma = NA_real_,
    gaussFitR2 = NA_real_,
    noNoiseObs = NA_integer_,
    noNoiseFraction = NA_real_,
    symmetry = NA_real_,
    signalBaselineRatio = NA_real_,
    groupOverlap = NA_real_,
    groupOverlapFrac = NA_real_,
    localMaxFlag = NA_real_,
    fromBlankSample = NA_integer_,
    label = NA_integer_
  )

  output <- list(
    peakgroups = peakgroups,
    samples = samples,
    peaks = peaks
  )

  check_defaults <- purrr::map(
    names(output),
    function(x) {
      tibble::tibble(
        table = x,
        variable = colnames(output[[x]]),
        current_class = unname(purrr::map_chr(output[[x]], class))
      )
    }
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::left_join(
      mzroll_schema_required_classes(),
      by = c("table", "variable")
    )

  class_mismatches <- check_defaults %>%
    dplyr::filter(current_class != required_class)
  if (nrow(class_mismatches) > 0) {
    stop(glue::glue("{nrow(class_mismatches)} defaults have the wrong class"))
  }

  unknown_defaults <- check_defaults[is.na(check_defaults$variable_is_required), ]
  if (nrow(unknown_defaults) > 0) {
    stop(glue::glue(
      "{nrow(unknown_defaults) defaults are not defined in
         clamr::mzroll_schema_required_classes"
    ))
  }

  return(output)
}


test_mzroll_db_con_schema_pks <- function(mzroll_db_con) {
  duplicated_peakgroups <- dplyr::tbl(
    mzroll_db_con,
    dbplyr::sql("SELECT peakgroups.groupId, COUNT(*) c
                 FROM peakgroups
                 GROUP BY peakgroups.groupId HAVING c > 1")
  ) %>%
    dplyr::collect()

  if (nrow(duplicated_peakgroups) != 0) {
    stop(nrow(duplicated_peakgroups), " groupIds were duplicated in \"peakgroups\"; each \"groupId\" must be unique")
  }

  duplicated_samples <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT samples.sampleId, COUNT(*) c
                                                              FROM samples
                                                              GROUP BY samples.sampleId HAVING c > 1")) %>%
    dplyr::collect()

  if (nrow(duplicated_samples) != 0) {
    stop(nrow(duplicated_samples), " sampleIds were duplicated in \"samples\"; each \"sampleId\" must be unique")
  }

  return(invisible(0))
}

#' Reduce mzroll samples
#'
#' Remove all sample-associated features (samples, scans, peaks) which are not from the samples in retained_sampleIds
#'
#' @inheritParams test_mzroll_db_con_schema
#' @param retained_sampleIds a vector of sampleId integers
#'
#' @return update sql database return 0 invisibly
#'
#' @export
reduce_mzroll_samples <- function(mzroll_db_con, retained_sampleIds) {
  stopifnot(class(retained_sampleIds) == "integer")

  present_tables <- DBI::dbListTables(mzroll_db_con)
  updated_tables <- c("samples", "peaks", "scans", "mz_update_key", "rt_update_key")

  for (a_table in updated_tables) {
    if (!(a_table %in% present_tables)) {
      warning(a_table, " is missing")
      next
    }

    sql_query <- paste0("SELECT * FROM ", a_table)

    updated_table <- dplyr::tbl(mzroll_db_con, dbplyr::sql(sql_query)) %>%
      dplyr::collect() %>%
      dplyr::filter(sampleId %in% retained_sampleIds)

    DBI::dbWriteTable(mzroll_db_con, a_table, updated_table, overwrite = TRUE)
  }

  # remove peakgroups which no longer contain peaks

  retained_groupIds <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT DISTINCT groupId FROM peaks")) %>%
    dplyr::collect()

  updated_peakgroups <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT * FROM peakgroups")) %>%
    dplyr::collect() %>%
    dplyr::semi_join(retained_groupIds, by = "groupId")

  DBI::dbWriteTable(mzroll_db_con, "peakgroups", updated_peakgroups, overwrite = TRUE)

  # shrink database to filtered data
  DBI::dbSendQuery(mzroll_db_con, "VACUUM")

  invisible(0)
}

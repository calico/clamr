#' Formula Monoisotopic Mass
#'
#' @param formulas a character vector of molecular formulas
#'
#' @return a tibble containing formula and exact_mass
#'
#' @examples
#' formulas <- c("C6H12O6", "C21H27N7O14P2", "C10H17N3O6S")
# " formula_monoisotopic_mass(formulas)
#'
#' @export
formula_monoisotopic_mass <- function(formulas) {
  stopifnot("character" %in% class(formulas))

  parsed_formulas <- parse_standard_formulas(unique(formulas))

  primary_isotopes <- clamr::isotope_summaries %>%
    dplyr::group_by(element) %>%
    dplyr::arrange(desc(isotope_freq)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  parsed_formulas %>%
    dplyr::left_join(primary_isotopes, by = "element") %>%
    dplyr::group_by(formula) %>%
    dplyr::summarize(exact_mass = sum(n * amu))
}

#' Parse formulas
#'
#' Separate chemical formulas into atom counts and return as a tibble.
#'
#' @param formulas a character vector of character formulas
parse_standard_formulas <- function(formulas) {
  standard_formulas <- grepl("^[A-Z][0-9A-Za-z]*$", formulas)
  if (sum(!standard_formulas) != 0) {
    warning(paste0(sum(!standard_formulas), " formulas are not standard and will be discarded (they do not match ^[0-9A-Za-z]+)"))
  }

  formulas <- formulas[standard_formulas]

  lapply(formulas, function(a_formula) {
    split_formula <- strsplit(a_formula, split = "")[[1]]
    running_formula <- NULL

    for (i in seq_along(split_formula)) {
      if (grepl("[A-Z]", split_formula[i])) {
        running_formula <- rbind(
          running_formula,
          data.frame(formula = a_formula, element = split_formula[i], n = "", stringsAsFactors = FALSE)
        )
      } else if (grepl("[a-z]", split_formula[i])) {
        running_formula[nrow(running_formula), "element"] <- paste0(running_formula[nrow(running_formula), "element"], split_formula[i])
      } else if (grepl("[0-9]", split_formula[i])) {
        running_formula[nrow(running_formula), "n"] <- paste0(running_formula[nrow(running_formula), "n"], split_formula[i])
      } else {
        warning(split_formula[i], " was not matched")
      }
    }
    running_formula$n <- as.numeric(sub("^$", 1, running_formula$n))

    running_formula
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::as_tibble()
}

#' Adduct Mass
#'
#' From a set of standard adduct formulas, e.g., [M+H]-, calculate the adducts mass and charge.
#'
#' @param adducts a character vector of adduct formulas
#'
#' @return a tibble containing one row per unique adduct
#'
#' @examples
#' adduct_mass(adducts = c("[M+H]+", "[M-NH4]2-", "[M]+"))
#'
#' @export
adduct_mass <- function(adducts) {
  checkmate::assertClass(adducts, "character")

  parsed_formulas <- parse_adducts(unique(adducts))

  primary_isotopes <- clamr::isotope_summaries %>%
    dplyr::group_by(element) %>%
    dplyr::arrange(desc(isotope_freq)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  parsed_formulas %>%
    dplyr::left_join(primary_isotopes, by = "element") %>%
    dplyr::group_by(adduct) %>%
    dplyr::summarize(
      adduct_mass = sum(n * amu),
      adduct_charge = -1 * n[element == "e"]
    )
}

parse_adducts <- function(adducts) {
  adduct_tbl <- tibble::tibble(adduct = unique(adducts))

  parsed_adducts <- adduct_tbl %>%
    # dplyr::filter(stringr::str_detect(adduct, adduct_pattern)) %>%
    dplyr::mutate(
      formula_change = purrr::map_chr(adduct, function(x) {
        stringr::str_match(x, "\\[M(.+)\\]")[2]
      }),
      charge = purrr::map_chr(adduct, function(x) {
        stringr::str_match(x, "\\](.+)$")[2]
      })
    ) %>%
    dplyr::mutate(
      formula_change_tbl = purrr::map(formula_change, parse_adduct_changes),
      formula_charge_tbl = purrr::map(charge, parse_charge_formulas)
    ) %>%
    dplyr::filter(!purrr::map_lgl(formula_charge_tbl, is.null)) %>%
    dplyr::select(adduct, formula_change_tbl, formula_charge_tbl) %>%
    tidyr::gather(tbl_type, tbl_value, -adduct) %>%
    tidyr::unnest(tbl_value) %>%
    dplyr::select(-tbl_type, -formula)

  return(parsed_adducts)
}

parse_adduct_changes <- function(formula) {
  formula_split <- split_formula(formula)

  if (is.null(formula_split)) {
    return(NULL)
  }

  formula_split_tbl <- tibble::tibble()
  for (formula_type in names(formula_split)) {
    partial_formula <- formula_split[[formula_type]]

    if (partial_formula != "") {
      parsed_formula <- suppressWarnings(parse_standard_formulas(partial_formula))
      if (nrow(parsed_formula) == 0) {
        parsed_formula <- suppressWarnings(parse_simple_formulas(partial_formula))

        if (nrow(parsed_formula) == 0) {
          warning(formula, " cound not be parsed, returning NULL")
          return(NULL)
        }
      }

      if (formula_type == "losses") {
        parsed_formula$n <- -1 * parsed_formula$n
      }

      formula_split_tbl <- dplyr::bind_rows(formula_split_tbl, parsed_formula)
    }
  }

  return(formula_split_tbl %>%
    dplyr::group_by(element) %>%
    dplyr::summarize(n = sum(n)))
}


parse_simple_formulas <- function(formulas) {
  formulas <- unique(formulas)
  standard_formulas <- grepl("^[0-9]*[A-Z][a-z]?$", formulas)
  if (sum(!standard_formulas) != 0) {
    warning(paste0(sum(!standard_formulas), " formulas are not standard and will be discarded (they do not match '^[0-9]*[A-Z][a-z]?')"))
  }

  tibble::tibble(formula = formulas[standard_formulas]) %>%
    dplyr::mutate(
      element = stringr::str_extract(formula, "[A-Z][a-z]?"),
      n = dplyr::case_when(
        stringr::str_detect(formula, "^[A-Z]") ~ 1L,
        TRUE ~ as.integer(stringr::str_extract(formula, "^[0-9]"))
      )
    )
}

parse_charge_formulas <- function(charge_strings) {
  stopifnot(class(charge_strings) %in% c(NULL, "character"))
  charge_strings <- unique(charge_strings[!(is.na(charge_strings) | charge_strings == "")])

  if (length(charge_strings) == 0) {
    return(NULL)
  }

  standard_charge_formulas <- grepl("^[0-9]?[+-]$", charge_strings)
  if (sum(!standard_charge_formulas) != 0) {
    warning(paste0(sum(!standard_charge_formulas), " formulas are not standard and will be discarded (they do not match '^[0-9]*[A-Z][a-z]?')"))
  }

  tibble::tibble(formula = charge_strings[standard_charge_formulas]) %>%
    dplyr::mutate(
      element = "e",
      n = dplyr::case_when(
        stringr::str_detect(formula, "^[0-9]") ~ as.integer(stringr::str_extract(formula, "^[0-9]")),
        TRUE ~ 1L
      )
    ) %>%
    dplyr::mutate(gain_multiplier = dplyr::case_when(
      stringr::str_detect(formula, "\\-") ~ 1,
      TRUE ~ -1
    )) %>%
    dplyr::mutate(n = n * gain_multiplier) %>%
    dplyr::select(-gain_multiplier)
}

#' Create chemical formula
#'
#' @param labeled_counts a tibble (or data.frame) containing two variables which will be used: label and count
#'
#' @return a length one character vector containing a chemical formula
create_isotopic_chemical_formula <- function(labeled_counts) {
  mass_gains <- labeled_counts %>%
    dplyr::filter(count > 0) %>%
    dplyr::mutate(count_label = ifelse(count == 1, label, paste0(label, "(", count, ")")))

  if (nrow(mass_gains) != 0) {
    gain_formula <- paste(mass_gains$count_label, collapse = " ")
  }

  mass_losses <- labeled_counts %>%
    dplyr::filter(count < 0) %>%
    dplyr::mutate(count_label = ifelse(count == -1, label, paste0(label, "(", -1 * count, ")")))

  if (nrow(mass_losses) != 0) {
    loss_formula <- paste(mass_losses$count_label, collapse = " ")
  }

  #' Explicitly check number of rows in mass_gains and mass_losses before
  #' returning.
  if (nrow(mass_gains) == 0 && nrow(mass_losses) == 0) {
    return("")
  } else if (nrow(mass_gains) != 0 && nrow(mass_losses) == 0) {
    return(gain_formula)
  } else if (nrow(mass_gains) == 0 && nrow(mass_losses) != 0) {
    return(paste("-", loss_formula))
  } else {
    return(paste(gain_formula, "-", loss_formula))
  }
}

#' Parse Isotopic Chemical Formula
#'
#' Read an internal, non-standard formula convention which allows for gains, losses and isotopes
#'
#' @param formula a chemical formula which may included losses as generated by \code{\link{create_isotopic_chemical_formula}}
#'
#' @rdname create_isotopic_chemical_formula
#'
#' @examples
#' parse_isotopic_chemical_formula("C(6) H(12) O(6)")
#' parse_isotopic_chemical_formula("13C C(5) H(12) O(6)")
#' parse_isotopic_chemical_formula("13C C(5) H(12) O(6) - H(2)")
#'
#' @export
parse_isotopic_chemical_formula <- function(formula) {
  formula_split <- split_formula(formula, collapse = " ")

  if (is.null(formula_split)) {
    return(NULL)
  }

  if (!is.null(formula_split$gains)) {
    gain_formula <- strsplit(formula_split$gains, split = " ")[[1]] %>%
      tibble::tibble(entry = ., is_positive = TRUE)
  }
  if (!is.null(formula_split$losses)) {
    loss_formula <- strsplit(formula_split$losses, split = " ")[[1]] %>%
      tibble::tibble(entry = ., is_positive = FALSE)
  }

  dplyr::bind_rows(gain_formula, loss_formula) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      element = regmatches(entry, regexpr("^[0-9A-Za-z]+", entry)),
      count = as.numeric(gsub("\\(|\\)", "", regmatches(entry, regexpr("\\([0-9]+\\)", entry)))) %>%
        ifelse(length(.) == 1, ., 1),
      count = ifelse(is_positive, count, count * -1)
    ) %>%
    dplyr::mutate(formula = formula) %>%
    dplyr::select(formula, element, count) %>%
    dplyr::group_by(formula, element) %>%
    dplyr::summarize(count = sum(count))
}

#' Split formula
#'
#' Split a formula into individual components
#'
#' @param formula a chemical formula which may included losses as generated by \code{\link{create_isotopic_chemical_formula}}
#' @param collapse delimiter used for collapsing individual formula into single string
#'
#' @export
split_formula <- function(formula, collapse = "") {
  if (is.na(formula)) {
    return(NULL)
  }

  gain_formula <- stringr::str_replace_all(formula, "\\- ?[0-9A-Za-z\\(\\) ]+", "") %>%
    stringr::str_replace("\\+", "") %>%
    stringr::str_trim() %>%
    stringr::str_replace_all(" ", collapse) %>%
    paste(collapse = collapse)

  # losses always start with a "-"
  loss_formula <- stringr::str_extract_all(formula, "\\- ?[0-9A-Za-z\\(\\) ]+")[[1]] %>%
    stringr::str_replace("-", "") %>%
    stringr::str_trim() %>%
    stringr::str_replace_all(" ", collapse) %>%
    paste(collapse = collapse)

  list(
    gains = gain_formula,
    losses = loss_formula
  )
}

#' Estimate chemical formula
#'
#' @param query_mass amu of an unknown specie
#' @param isotope_summary isotopes that the unknown specie could be constructed with:
#' \itemize{
#'   \item{label}
#'   \item{amu}
#'   \item{count_lb (optional)}
#'   \item{count_ub (optional)}
#'   \item{prior_weight (optional)}}
#' If penalization = "prior", a variable named "prior" must be specified (larger weights are bad)".
#' @param penalization the method of penalizing the optimization
#' \itemize{
#'  \item{none: minimize the absolute difference between predicted mass and "query_mass"}
#'  \item{count: none + minimize sum of absolute atom counts}
#'  \item{prior: none + minimize a weighted sum of absolute atom counts}}
#' @param penalization_multiplier multiplier for the penalization functions within the cost function
#' @param timeout seconds until a timeout occurs
#'
#' @return a list containing a summary of fit (absolute mass difference and penalization cost (if penalization != "none")) and isotope_summary with an added variable
#'
#' @examples
#' library(dplyr)
#'
#' probable_isotopes <- clamr::isotope_summaries %>%
#'   dplyr::inner_join(
#'     clamr::elemental_frequency_summary %>%
#'       dplyr::filter(human_body > 0.0001) %>%
#'       dplyr::select(element, human_body),
#'     by = "element"
#'   ) %>%
#'   dplyr::mutate(natural_freq = isotope_freq * human_body) %>%
#'   dplyr::filter(natural_freq > 0.0001) %>%
#'   dplyr::mutate(
#'     count_lb = ifelse(element == "H", -8, -2),
#'     count_lb = ifelse(element %in% c("C", "O"), -4, count_lb)
#'   ) %>%
#'   dplyr::mutate(count_ub = -1 * count_lb) %>%
#'   dplyr::rowwise() %>%
#'   dplyr::mutate(label = paste0(element, "-", strsplit(as.character(amu), split = "\\.")[[1]][1])) %>%
#'   dplyr::ungroup() %>%
#'   dplyr::mutate(prior_weight = pmin(1 / 100 / natural_freq, 1))
#'
#' query_mass <- -1.007825 + 13.003355 + 2 * 12
#' estimate_chemical_formula(query_mass, probable_isotopes, penalization = "count")
#'
#' @export
estimate_chemical_formula <- function(query_mass, isotope_summary, penalization = "none", penalization_multiplier = 0.001, timeout = 10) {

  # check query_mass

  stopifnot(length(query_mass) == 1, class(query_mass) == "numeric")

  # check formatting of isotope_summary

  stopifnot("data.frame" %in% class(isotope_summary), all(c("label", "amu") %in% colnames(isotope_summary)))

  reserved_variables <- "count"
  if (any(reserved_variables %in% colnames(isotope_summary))) {
    stop(paste(reserved_variables, collapse = " & "), " are reserved variables in isotope_summary")
  }

  # check penalization

  stopifnot(length(penalization) == 1, class(penalization) == "character")
  penalization_options <- c("none", "count", "prior")
  if (!(penalization %in% penalization_options)) {
    stop(paste0('"penalization" has an invalid value; valid values are: ', paste(penalization_options, collapse = ", ")))
  }

  if (penalization == "prior" & !("prior_weight" %in% colnames(isotope_summary))) {
    stop('No prior weight found: if penalizaation equals "prior", a variable "prior_weight" must be provided in "isotope_summary"')
  }

  # check penalizaiton_multiplier

  stopifnot(length(penalization_multiplier) == 1, class(penalization_multiplier) == "numeric", penalization_multiplier >= 0)

  # check for whether isotope counts are bounded
  if (all(c("count_lb", "count_ub") %in% colnames(isotope_summary))) {
    counts_are_bounded <- TRUE
  } else {
    counts_are_bounded <- FALSE
    message('counts will not be bounded because "count_lb" and "count_ub" were not present in isotope_summary')
  }

  n_species <- nrow(isotope_summary)
  specie_masses <- isotope_summary$amu

  # estimate chemical formula as a mixed integer linear programming (MILP) optimization problem

  if (penalization == "none") {

    # minimize the absolute deviation between query mass and predicted chemical formula mass

    lp_model <- lpSolveAPI::make.lp(nrow = 0, ncol = n_species + 3)

    # separate objective into a negative and positive component with delta variables indicating the positive or negative departure
    lpSolveAPI::add.constraint(lp_model, c(isotope_summary$amu, -1 * query_mass, 1, 0), ">=", 0)
    lpSolveAPI::add.constraint(lp_model, c(isotope_summary$amu, -1 * query_mass, 0, 1), "<=", 0)

    lpSolveAPI::set.objfn(lp_model, c(rep(0, n_species), 0, 1, -1))
    # all counts of species are integers
    lpSolveAPI::set.type(lp_model, columns = 1:n_species, type = "integer")
    # constraint on measured mz
    lpSolveAPI::set.bounds(lp_model, lower = 1, upper = 1, columns = n_species + 1)
    # add absolute value constraint
    lpSolveAPI::set.bounds(lp_model, lower = 0, upper = Inf, columns = n_species + 2)
    lpSolveAPI::set.bounds(lp_model, lower = -Inf, upper = 0, columns = n_species + 3)
    # optinally constrain counts
    if (counts_are_bounded) {
      lpSolveAPI::set.bounds(lp_model, lower = isotope_summary$count_lb, upper = isotope_summary$count_ub, columns = 1:n_species)
    }

    lpSolveAPI::lp.control(lp_model, timeout = timeout, epsb = 1e-20, epsd = 1e-20, epsint = 1e-10) # need a low integer tolerance
    status <- solve(lp_model)

    objective <- tibble::tibble(
      query_mass = query_mass,
      objective = lpSolveAPI::get.objective(lp_model),
      status = status,
      absolute_mass_diff = lpSolveAPI::get.objective(lp_model)
    )
  } else {
    # minimize the absolute deviation between query mass and predicted chemcial formula mass + regularization based on fxn(#/types of atoms)

    lp_model <- lpSolveAPI::make.lp(nrow = 0, ncol = n_species * 3 + 4)
    zero_vector <- rep(0, n_species * 3 + 4)

    # separate objective into a negative and positive component with delta variables indicating the positive or negative departure
    lpSolveAPI::add.constraint(lp_model, c(isotope_summary$amu, rep(0, n_species * 2), -1 * query_mass, 1, 0, 0), ">=", 0)
    lpSolveAPI::add.constraint(lp_model, c(isotope_summary$amu, rep(0, n_species * 2), -1 * query_mass, 0, 1, 0), "<=", 0)

    lpSolveAPI::set.objfn(lp_model, c(rep(0, n_species * 3), 0, 1, -1, penalization_multiplier))
    # all counts of species are integers
    lpSolveAPI::set.type(lp_model, columns = 1:n_species, type = "integer")
    # constraint on measured mz
    lpSolveAPI::set.bounds(lp_model, lower = 1, upper = 1, columns = n_species * 3 + 1)
    # add absolute value constraint
    lpSolveAPI::set.bounds(lp_model, lower = 0, upper = Inf, columns = n_species * 3 + 2)
    lpSolveAPI::set.bounds(lp_model, lower = -Inf, upper = 0, columns = n_species * 3 + 3)
    # optinally constrain counts
    if (counts_are_bounded) {
      lpSolveAPI::set.bounds(lp_model, lower = isotope_summary$count_lb, upper = isotope_summary$count_ub, columns = 1:n_species)
    }

    # add the absolute value of counts (as a positive components plus -1*a negative component)
    # both components are positive -- one component will be nonzero and a summation variable will equal this value
    for (i in 1:n_species) {
      # negative component
      assign_vec <- zero_vector
      assign_vec[c(i, n_species + i)] <- c(1, -1)
      lpSolveAPI::add.constraint(lp_model, assign_vec, "<=", 0)
      # positive component
      assign_vec <- zero_vector
      assign_vec[c(i, n_species * 2 + i)] <- c(1, 1)
      lpSolveAPI::add.constraint(lp_model, assign_vec, ">=", 0)
      # sum of both components
      assign_vec <- zero_vector
      assign_vec[c(n_species + i, n_species * 2 + i)] <- c(1, 1)
      lpSolveAPI::add.constraint(lp_model, assign_vec, ">=", 0)
    }
    # sum or prior on atom counts
    if (penalization == "count") {
      assign_vec <- zero_vector
      assign_vec[(n_species + 1):(n_species * 3)] <- 1
      assign_vec[length(assign_vec)] <- -1
    } else if (penalization == "prior") {
      assign_vec <- zero_vector
      assign_vec[(n_species + 1):(n_species * 3)] <- rep(isotope_summary$prior_weight, times = 2)
      assign_vec[length(assign_vec)] <- -1
    } else {
      stop("penalization type not recognized")
    }
    lpSolveAPI::add.constraint(lp_model, assign_vec, "=", 0)

    lpSolveAPI::lp.control(lp_model, timeout = timeout, epsb = 1e-20, epsd = 1e-20, epsint = 1e-10)
    status <- solve(lp_model)

    costs <- lpSolveAPI::get.variables(lp_model)[n_species * 3 + 2:4]
    objective <- tibble::tibble(
      query_mass = query_mass,
      objective = lpSolveAPI::get.objective(lp_model),
      status = status,
      absolute_mass_diff = sum(costs[1:2]),
      penalization = costs[3] * penalization_multiplier
    )
  }

  # output objective (split up if there are multiple objectives) and associated formula

  formula_table <- isotope_summary %>%
    dplyr::mutate(
      count = lpSolveAPI::get.variables(lp_model)[1:n_species],
      query_mass = query_mass
    ) %>%
    dplyr::select_(.dots = c("query_mass", "label", "count", intersect("prior_weight", colnames(isotope_summary)))) %>%
    dplyr::filter(count != 0)

  objective$formula <- create_isotopic_chemical_formula(formula_table)

  list(
    objective = objective,
    formula = formula_table
  )
}

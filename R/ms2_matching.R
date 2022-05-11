#' Score Fragmentation Similarity
#'
#' Score a single possible ion match based on the agreement of an experimental and standard spectra.
#'
#' @export
#'
#' @param an_expFrag a tibble of experimental fragmentation spectra data for a single ion/energy.
#' @param a_stdFrag a tibble of standard fragmentation spectra data for a single ion/energy.
#' @param MS2tol the instrument mass accuracy tolerance for the purpose of fragment grouping - \code{\link[clamr]{format_mass_accuracy_input}}, normally passed from a \code{\link[clamr]{build_clamr_config}} object.
#' @param frag_similarity_methods a character vector of methods to use for fragmentation scoring:
#' \itemize{
#'   \item{tic: fraction of experimental signal explained by standard fragments.}
#'   \item{fisher: fisher exact test to look at overlap between experimental and standard fragments.}
#'   \item{spearman: spearman correlation between experimental and standard fragments.}
#'   \item{cosine: the cosine similarity between mz-matched ions of the standard and experimental spectra.}
#'   \item{robust cosine: calculate the cosine similarity when separately dropping out each fragment - then take the worst cosine similarity among them. Identifies cases where a match score is driven by a single fragment.}
#'  }
#' @param remove_precursor Should the precursor mass be removed before matching spectra?
#' @param minimum_ic_fraction A value between 0 and 1 indicating for a given peak the fraction of total signal for it to be used for matching.
#' @param print_plot Should a summary plot be printed?
#'
#' @return a tibble containing frag_similarity_method and method_score for this match.
score_fragmentation_similarity <- function(an_expFrag, a_stdFrag, MS2tol, frag_similarity_methods = "cosine", remove_precursor = FALSE, minimum_ic_fraction = 0, print_plot = FALSE) {
  available_similarity_methods <- tibble::frame_data(
    ~method_flag, ~method_fxn,
    "tic", "score_fragmentation_similarity_tic_percent",
    "fisher", "score_fragmentation_similarity_fisher_exact",
    "spearman", "score_fragmentation_similarity_spearman",
    "cosine", "score_fragmentation_similarity_cosine",
    "robust_cosine", "score_fragmentation_similarity_robust_cosine"
  )

  stopifnot(class(frag_similarity_methods) == "character")
  if (any(frag_similarity_methods == "all")) {
    frag_similarity_methods <- available_similarity_methods$method_flag
  }

  undefined_methods <- setdiff(frag_similarity_methods, available_similarity_methods$method_flag)
  if (length(undefined_methods) != 0) {
    warning(paste(undefined_methods, collapse = " & "), " are not defined frag_similarity_methods")
  }
  usable_methods <- intersect(available_similarity_methods$method_flag, frag_similarity_methods)
  if (length(usable_methods) == 0) {
    stop("no valid methods were supplied in frag_similarity_methods:, valid options are ", paste(available_similarity_methods$method_flag, collapse = " & "))
  }
  stopifnot(length(remove_precursor) == 1, class(remove_precursor) == "logical", remove_precursor %in% c(TRUE, FALSE))
  stopifnot(length(minimum_ic_fraction) == 1, class(minimum_ic_fraction) == "numeric", minimum_ic_fraction >= 0, minimum_ic_fraction < 1)
  stopifnot(length(print_plot) == 1, class(print_plot) == "logical", print_plot %in% c(TRUE, FALSE))

  # apply filters

  if (remove_precursor) {
    if ("is_precursor" %in% colnames(an_expFrag) && "is_precursor" %in% colnames(a_stdFrag)) {
      an_expFrag <- an_expFrag %>% dplyr::filter(!is_precursor)
      a_stdFrag <- a_stdFrag %>% dplyr::filter(!is_precursor)
    } else {
      warning("\"remove_precursor\" is TRUE, but \"is_precursor\" was not found in an_expFrag or a_stdFrag,\nprecursors will not be removed.\n")
    }
  }

  # filter low intensity fragments

  if (minimum_ic_fraction != 0) {
    an_expFrag <- spectra_filter_fragments(an_expFrag, minimum_ic_fraction)
    a_stdFrag <- spectra_filter_fragments(a_stdFrag, minimum_ic_fraction)

    if (is.null(an_expFrag) || is.null(a_stdFrag)) {
      return(NULL)
    }
  }

  # print a summary plot if desired
  if (print_plot) {
    print(plot_score_fragmentation_similarity(an_expFrag, a_stdFrag, MS2tol))
  }

  # score matches

  lapply(usable_methods, function(a_method_flag) {
    method_fxn <- available_similarity_methods$method_fxn[available_similarity_methods$method_flag == a_method_flag]
    # match optional method args to dots
    method_score <- do.call(method_fxn, list(
      an_expFrag = an_expFrag,
      a_stdFrag = a_stdFrag,
      MS2tol = MS2tol
    ))

    tibble::frame_data(
      ~frag_similarity_method, ~method_score,
      a_method_flag, method_score
    )
  }) %>%
    dplyr::bind_rows()
}

#' Spectra Filter Fragments
#'
#' @param fragData a tibble of fragmentation data such as those used in \code{\link{score_fragmentation_similarity}}.
#' @param minimum_ic_fraction minimum fraction of total IC for a fragment to be retained.
#'
#' @return \code{fragData} with low intensity fragments filtered.
#'
#' @export
spectra_filter_fragments <- function(fragData, minimum_ic_fraction = 0.001) {
  stopifnot(class(minimum_ic_fraction) == "numeric", length(minimum_ic_fraction) == 1, minimum_ic_fraction >= 0, minimum_ic_fraction < 1)

  if (minimum_ic_fraction > 0.1) {
    warning("minimum_ic_fraction is set at a high level: ", minimum_ic_fraction, " - most fragmentation peaks may be filtered\n")
  }

  fragData <- fragData %>%
    dplyr::filter(ic >= sum(ic) * minimum_ic_fraction)

  if (nrow(fragData) == 0) {
    warning("all spectra fragments filtered - \"minimum_ic_fraction\" may be set too high in \"spectra_filter_fragments\"\n")
    return(NULL)
  }

  fragData
}

plot_score_fragmentation_similarity <- function(an_expFrag, a_stdFrag, MS2tol) {
  binned_spectra <- combine_exp_std_fragmentations(an_expFrag, a_stdFrag, MS2tol)

  mz_labels <- binned_spectra %>%
    dplyr::group_by(mz_bin) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(mzlabel = as.character(round(mz, digits = 2))) %>%
    dplyr::ungroup() %>%
    {
      structure(.$mz_bin, names = .$mzlabel)
    }

  binned_spectra %>%
    dplyr::group_by(type) %>%
    dplyr::mutate(ic_frac = ic / sum(ic)) %>%
    ggplot(aes(x = mz_bin, y = ic_frac, fill = type)) +
    geom_bar(stat = "identity", position = "stack", width = 0.35) +
    scale_x_discrete("Fragment M/Z", limits = mz_labels, labels = names(mz_labels)) +
    scale_y_continuous("Fraction of total signal", labels = scales::percent) +
    scale_fill_brewer("Spectra type", palette = "Set2") +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(angle = 90, hjust = 1)
    )
}

score_fragmentation_similarity_tic_percent <- function(an_expFrag, a_stdFrag, MS2tol) {
  experimental_ic_total <- sum(an_expFrag$ic)

  experimental_ic_explained <- combine_exp_std_fragmentations(an_expFrag, a_stdFrag, MS2tol) %>%
    dplyr::group_by(mz_bin) %>%
    dplyr::filter(any("standard" %in% type)) %>%
    dplyr::filter(type == "experimental") %>%
    dplyr::ungroup() %>%
    {
      sum(.$ic)
    }

  experimental_ic_explained / experimental_ic_total
}

score_fragmentation_similarity_fisher_exact <- function(an_expFrag, a_stdFrag, MS2tol) {
  cont_table <- combine_exp_std_fragmentations(an_expFrag, a_stdFrag, MS2tol) %>%
    dplyr::distinct(mz_bin, type) %>%
    dplyr::mutate(present = 1) %>%
    tidyr::spread(type, present, fill = 0) %>%
    dplyr::mutate(
      experimental = factor(experimental, levels = c(0, 1)),
      standard = factor(standard, levels = c(0, 1))
    ) %>%
    {
      table(.$experimental, .$standard)
    }

  fragment_padding <- pmax(1000, sum(cont_table) * 2)

  cont_table[1, 1] <- fragment_padding - sum(cont_table)
  -1 * log10(stats::fisher.test(cont_table)$p.value)
}

score_fragmentation_similarity_spearman <- function(an_expFrag, a_stdFrag, MS2tol) {
  combine_exp_std_fragmentations(an_expFrag, a_stdFrag, MS2tol) %>%
    dplyr::group_by(mz_bin, type) %>%
    dplyr::summarize(ic = sum(ic)) %>%
    tidyr::spread(type, ic, fill = 0) %>%
    dplyr::ungroup() %>%
    dplyr::select(-mz_bin) %>%
    as.matrix() %>%
    {
      stats::cor(.[, 1], .[, 2], method = "spearman")
    }
}

score_fragmentation_similarity_cosine <- function(an_expFrag, a_stdFrag, MS2tol) {
  combine_exp_std_fragmentations(an_expFrag, a_stdFrag, MS2tol) %>%
    dplyr::group_by(mz_bin, type) %>%
    dplyr::summarize(ic = sum(ic)) %>%
    tidyr::spread(type, ic, fill = 0) %>%
    dplyr::ungroup() %>%
    dplyr::select(-mz_bin) %>%
    as.matrix() %>%
    {
      sum(.[, "experimental"] * .[, "standard"]) / sqrt(sum(.[, "experimental"]^2) * sum(.[, "standard"]^2))
    }
}

score_fragmentation_similarity_robust_cosine <- function(an_expFrag, a_stdFrag, MS2tol) {
  spread_ions <- combine_exp_std_fragmentations(an_expFrag, a_stdFrag, MS2tol) %>%
    dplyr::group_by(mz_bin, type) %>%
    dplyr::summarize(ic = sum(ic)) %>%
    tidyr::spread(type, ic, fill = 0)

  # separately drop each mz and calculate cosine similarity

  robust_cosine_similarities <- expand.grid(held_out = spread_ions$mz_bin, mz_bin = spread_ions$mz_bin, stringsAsFactors = FALSE) %>%
    dplyr::filter(held_out != mz_bin) %>%
    dplyr::left_join(spread_ions, by = "mz_bin") %>%
    dplyr::group_by(held_out) %>%
    dplyr::summarize(robust_cosine = sum(experimental * standard) / sqrt(sum(experimental^2) * sum(standard^2))) %>%
    dplyr::arrange(robust_cosine)

  robust_cosine_similarities$robust_cosine[1]
}

#' @inheritParams score_fragmentation_similarity
combine_exp_std_fragmentations <- function(an_expFrag, a_stdFrag, MS2tol) {
  clamr::sum_spectra(a_stdFrag, MS2tol) %>%
    dplyr::select(mz, ic) %>%
    dplyr::mutate(type = "standard") %>%
    dplyr::bind_rows(
      an_expFrag %>%
        dplyr::mutate(type = "experimental")
    ) %>%
    dplyr::arrange(mz) %>%
    dplyr::mutate(mz_bin = clamr::find_mz_jumps(.$mz,
      MS2tol$tol,
      MS2tol$absolute_or_relative,
      collapse = FALSE
    ))
}

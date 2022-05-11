stratify_signal_frac <- function(group_data) {
  defined_group_data <- group_data %>%
    dplyr::filter(!is.na(signal_frac_minor))

  if (length(unique(defined_group_data$stratify_category)) == 1) {
    return(tibble::tibble(term = "stratify_category"))
  }

  defined_group_data %>%
    stats::aov(formula = signal_frac_minor ~ stratify_category) %>%
    broom::tidy() %>%
    dplyr::filter(term == "stratify_category")
}

cosine_similarity <- function(x) {
  x %*% t(x) / (sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
}

#' Separate by local minima
#'
#' Determine the density of x (using the density function from stats) and then bin all points x using local minima as boundaries
#'
#' @param x a numeric vector
#' @param bw the sd of the KDE smoother
#' @param print_plot TRUE/FALSE print a plot showing
#'
#' @return a vector of integer group labels
find_density_minima <- function(x, bw = NULL, print_plot = FALSE) {
  stopifnot(class(bw) %in% c("NULL", "numeric", "integer"))
  stopifnot(class(print_plot) == "logical", length(print_plot) == 1, print_plot %in% c(TRUE, FALSE))
  if (class(bw) == "NULL") {
    bw <- diff(range(x)) / 100
  } else {
    stopifnot(length(bw) == 1, bw > 0)
  }

  nbins <- diff(range(x)) / bw * 10

  if (nbins < 3) {
    # never split if there is very little variation in x
    return(rep(0, times = length(x)))
  }
  if (nbins > 2^20) {
    stop("the range of x is too great relative to bw; split x using other methods first")
  }

  stopifnot(class(print_plot) == "logical", length(print_plot) == 1, print_plot %in% c(TRUE, FALSE))

  x_density <- stats::density(x, bw = bw, n = nbins)
  which_minima <- (x_density$y[1:(nbins - 2)] > x_density$y[2:(nbins - 1)]) & (x_density$y[3:nbins] > x_density$y[2:(nbins - 1)])
  values_minima <- c(0, x_density$x[2:(nbins - 1)][which_minima])

  bin_assignment <- sapply(x, function(x) {
    sum(x > values_minima)
  })

  if (print_plot) {
    print(qplot(x = x, fill = as.factor(bin_assignment %% 2), bins = 1000) +
      theme_bw() + theme(text = element_text(size = 20)) +
      scale_fill_brewer(palette = "Set1", guide = FALSE))
  }

  bin_assignment
}

#' Merge Peakgroups
#'
#' Overwrites peakgroups and peaks tables in mzroll_db_con after merging
#' peakgroup pairs specified in split_peaks.
#'
#' @inheritParams test_mzroll_db_con_schema
#' @param split_peaks output of clamr::identify_split_peaks
#'
#' @return 0 invisibly (run for side effects)
#'
#' @export
merge_peakgroups <- function(mzroll_db_con, split_peaks) {

  # find all sub-graphs when treating split_peaks as an edge list
  # to collapse cases of A:B, B:C, C:D into a single new groupId
  group_updates <- igraph::graph_from_data_frame(split_peaks, directed = FALSE) %>%
    {
      tibble::tibble(
        groupId_old = as.integer(igraph::V(.)$name),
        group = igraph::clusters(.)$membership
      )
    } %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(groupId_new = min(groupId_old)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-group)

  peakset <- clamr::extract_peakset(mzroll_db_con) %>%
    peakset_format_for_mzroll()

  untouched_peaks <- peakset$peaks %>%
    dplyr::anti_join(group_updates, by = c("groupId" = "groupId_old"))
  updated_peaks <- peakset$peaks %>%
    dplyr::inner_join(group_updates, by = c("groupId" = "groupId_old")) %>%
    # for each group and sample choose the highest intensity peak
    dplyr::group_by(groupId_new, sampleId) %>%
    dplyr::arrange(desc(peakAreaTop)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-groupId) %>%
    dplyr::rename(groupId = groupId_new)
  peak_overwrite <- dplyr::bind_rows(untouched_peaks, updated_peaks)

  untouched_peakgroups <- peakset$peakgroups %>%
    dplyr::anti_join(group_updates, by = c("groupId" = "groupId_old"))
  updated_peakgroups <- peakset$peakgroups %>%
    dplyr::semi_join(group_updates, by = c("groupId" = "groupId_new"))
  peakgroup_overwrite <- dplyr::bind_rows(untouched_peakgroups, updated_peakgroups)

  message(glue::glue("{nrow(group_updates)} split peakgroups consolidated into {nrow(updated_peakgroups)} peakgroups"))

  clamshell::sqlite_write(mzroll_db_con, "peaks", peak_overwrite, overwrite = TRUE)
  clamshell::sqlite_write(mzroll_db_con, "peakgroups", peakgroup_overwrite, overwrite = TRUE)

  invisible(0)
}

#' Group Similar MS2
#'
#' Wraps \code{\link{find_consistent_ms2_fingerprints}} to find consistent ms2 groups with high information
#'
#' @inheritParams find_consistent_ms2_fingerprints
#' @inheritParams test_clamr_config
#' @inheritParams find_consistent_ms2_fingerprints
#'
#' @examples
#' clamr_config <- build_clamr_config(list(MS2tol = "20ppm"))
#' grouped_ms2s <- group_similar_MS2(mz_set = clamr::mz_set_example, clamr_config, cosine_cutoff = 0.95)
#'
#' dplyr::count(grouped_ms2s$ms2_groups, ms2_group)
#'
#' @export
group_similar_MS2 <- function(mz_set, clamr_config, cosine_cutoff = 0.95) {
  # cluster MS2 spectra with similar precursor mz

  ms2_groups <- find_consistent_ms2_fingerprints(mz_set,
    MS2tol = clamr_config$MS2tol,
    cosine_cutoff = cosine_cutoff
  )

  if (class(ms2_groups) == "NULL") {
    # no consistent fragments in spectra
    return(NULL)
  }

  list(
    ms2_groups = mz_set %>%
      dplyr::mutate(peak_num = 1:dplyr::n()) %>%
      dplyr::left_join(ms2_groups$ms2_group_assignments, by = "peak_num"),
    ms2_footprints = ms2_groups$ms2_footprints
  )
}

#' Find Consistent MS2 Fingerprints
#'
#' Compare all MS2 fingerprints within the group to separate distinct groups based on MS2 profiles.
#'
#' @param mz_set a set of peaks with a consistent parent mass m/z which may contain multiple species with distinct RT / MS2 fingerprints.
#' @param MS2tol mass tolerance for grouping ms2 m/z values: see \code{\link{build_clamr_config}}.
#' @param cosine_cutoff Minimum cosine similarity between a pair of MS2 fragment profiles to group them into a common cluster.
#' @param create_plots TRUE/FALSE: should a summary plot be printed.
#' @param rt_range optional if create_plots is FALSE, if create_plots is TRUE, length two numeric vector specifying the lower and upper limits of retention time for an experiment.
#'
#' @return assignments of parent peaks to distinct ms2 peak groups
#'
#' @examples
#' MS2tol <- build_clamr_config(list(MS2tol = "20ppm"))$MS2tol
#' find_consistent_ms2_fingerprints(
#'   mz_set = clamr::mz_set_example,
#'   MS2tol = MS2tol, create_plots = TRUE, rt_range = c(0, 22)
#' )
#'
#' @export
find_consistent_ms2_fingerprints <- function(mz_set, MS2tol, cosine_cutoff = 0.95, create_plots = FALSE, rt_range = NULL) {
  stopifnot(class(cosine_cutoff) == "numeric", length(cosine_cutoff) == 1, all(cosine_cutoff >= 0 & cosine_cutoff <= 1))
  stopifnot(class(create_plots) == "logical", create_plots %in% c(TRUE, FALSE))
  stopifnot(!(is.null(rt_range) & create_plots == TRUE))
  stopifnot(class(rt_range) %in% c("NULL", "numeric"))
  if (class(rt_range) == "numeric") {
    stopifnot(length(rt_range) == 2)
  }

  mz_set <- mz_set %>%
    dplyr::mutate(peak_num = 1:dplyr::n())

  # unnest all MS2 fragments and group fragments based on MS2 tolerance
  grouped_ms2_data <- extract_and_group_fragments(mz_set, MS2tol)

  # filter low information fragments
  filtered_grouped_ms2_data <- filter_low_information_fragments(grouped_ms2_data, filter_method = "signal_pool")

  # filter low complexity spectra
  filtered_grouped_ms2_data <- filter_low_complexity_spectra(filtered_grouped_ms2_data, min_signal_fragments = 5, signal_fraction_co = 0.01)

  if (length(unique(filtered_grouped_ms2_data$peak_num)) <= 1) {
    # if there are no consistent fragments or less than 1 parent peak with major fragments then return NULL
    return(NULL)
  }

  if (create_plots) {
    if (is.null(rt_range)) {
      rt_range <- range(mz_set$rt)
    }
    # optionally display MS2 fingerprints for each precursor peak, along with the RT and IC of the precursor for reference
    plot_ms1_group(mz_set, grouped_ms2_data,
      rt_range = rt_range,
      plot_title = paste0("Parent M/Z:", as.character(round(stats::median(mz_set$precursorMz), 4)))
    )
  }

  grouped_ms2_parent_similarity <- grouped_ms2_data %>%
    reshape2::acast(peak_num ~ mz_set, value.var = "ms2_ic", fill = 0) %>%
    cosine_similarity()
  cosine_peak_order <- as.character(fastcluster::hclust(stats::as.dist(grouped_ms2_parent_similarity), method = "single")$order)

  # group parent peaks based on ms2 similarity assuming a block-diagonal structure
  # by assuming a block-diagonal structure, I can look at the cor between i and i-1 to assesss correlation,
  # rather than using the average of a full block

  sample_adjacencies <- grouped_ms2_parent_similarity %>%
    # converting from matrix back to tidy data
    as.data.frame() %>%
    dplyr::mutate(row = rownames(.)) %>%
    tidyr::gather("col", "cos", -row) %>%
    # hierarchical clustering using cosine distance
    dplyr::mutate(
      row = factor(row, levels = cosine_peak_order),
      col = factor(col, levels = cosine_peak_order)
    ) %>%
    # coerce factors back to numerics to compare adjacent values that were adjacent in hclust
    dplyr::filter(as.numeric(row) == as.numeric(col) + 1) %>%
    # identify block-diagonal regions of similar MS2 spectra
    dplyr::mutate(block_diag = ifelse(cos > cosine_cutoff, TRUE, FALSE))

  # determine ms2 group membership from breaks in the block diagonal structure

  parent_asignments <- (sample_adjacencies %>%
    dplyr::filter(block_diag) %>%
    # convert to network
    igraph::graph_from_data_frame() %>%
    igraph::clusters())$membership

  parent_singletons <- setdiff(as.character(mz_set$peak_num), names(parent_asignments))

  if (length(parent_singletons) != 0) {
    parent_singletons <- tibble::tibble(peak_num = as.integer(parent_singletons)) %>%
      dplyr::mutate(ms2_group = 1:dplyr::n() + max(c(suppressWarnings(max(parent_asignments)), 0)))
  } else {
    parent_singletons <- NULL
  }

  if (length(parent_asignments) != 0) {
    parent_asignments <- tibble::tibble(peak_num = as.integer(names(parent_asignments)), ms2_group = unname(parent_asignments))
  } else {
    parent_asignments <- NULL
  }

  final_group_assignments <- dplyr::bind_rows(parent_asignments, parent_singletons)

  peak_profiles <- grouped_ms2_data %>%
    dplyr::left_join(final_group_assignments, by = "peak_num") %>%
    dplyr::group_by(ms2_group, ms2_mz) %>%
    dplyr::summarize(n_parents = dplyr::n(), mean_ms2_ic_frac = mean(ms2_ic_frac))

  list(
    ms2_group_assignments = final_group_assignments,
    ms2_footprints = peak_profiles
  )
}

#' Filter Low Information Fragments
#'
#' Filter to high information fragments based on one of several filter_methods
#'
#' @param grouped_ms2_data output of \code{\link{extract_and_group_fragments}}
#' @param filter_method method for filtering low-information fragments:
#' \describe{
#'   \item{top5union (default)}{the union of the 5 most abundant fragments (w/ intensity >1\% of total signal)}
#'   \item{frequency}{take fragments which are present in more than 10 spectra or >= 20\% of spectra or that are the most abundant ion of any spectra.}
#'   \item{signal_pool}{take fragments which represents more than \code{signal_fraction_co} fraction of the total signal across all peak spectra.}
#' }
#' @param signal_fraction_co fraction of total signal that a fragment must represent to be retained when using "signal_pool".
#'
#' @return grouped_ms2_data with low-information fragments removed.
filter_low_information_fragments <- function(grouped_ms2_data, filter_method = "top5union", signal_fraction_co = 0.001) {
  stopifnot("data.frame" %in% class(grouped_ms2_data))
  stopifnot(
    class(filter_method) == "character", length(filter_method) == 1,
    filter_method %in% c("top5union", "frequency", "signal_pool")
  )
  stopifnot(class(signal_fraction_co) == "numeric", length(signal_fraction_co) == 1, signal_fraction_co > 0, signal_fraction_co < 1)

  if (filter_method == "top5union") {
    # take the union of the top 5 most abundant ions in each spectra
    most_abundant_ms2_set <- grouped_ms2_data %>%
      dplyr::group_by(peak_num) %>%
      dplyr::filter(ms2_ic_frac > 0.01) %>%
      dplyr::arrange(desc(ms2_ic_frac)) %>%
      dplyr::slice(1:5) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(mz_set)

    grouped_ms2_data <- grouped_ms2_data %>%
      dplyr::group_by(mz_set) %>%
      dplyr::filter(mz_set %in% most_abundant_ms2_set$mz_set)
  } else if (filter_method == "frequency") {
    # take fragments which are present in more than 10 spectra or >= 20% of spectra or that are the most abundant ion of any spectra
    most_abundant_ms2_set <- grouped_ms2_data %>%
      dplyr::group_by(peak_num) %>%
      dplyr::filter(ms2_ic == max(ms2_ic)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(mz_set)

    # filter rare peaks
    grouped_ms2_data <- grouped_ms2_data %>%
      dplyr::group_by(mz_set) %>%
      # filter low signal MS2 fragments which are rare
      dplyr::filter(dplyr::n() > 10 | dplyr::n() / nrow(.) >= 0.2 | ms2_mz %in% most_abundant_ms2_set$ms2_mz)
  } else if (filter_method == "signal_pool") {
    # take fragments which represents more than 1% of the total signal across all peak spectra
    summed_signal_fractions <- grouped_ms2_data %>%
      dplyr::group_by(ms2_mz) %>%
      dplyr::summarize(signal_fraction_sum = sum(ms2_ic_frac)) %>%
      dplyr::filter(signal_fraction_sum > signal_fraction_co * length(unique(grouped_ms2_data$peak_num)))

    grouped_ms2_data <- grouped_ms2_data %>%
      dplyr::filter(ms2_mz %in% summed_signal_fractions$ms2_mz)
  }

  return(grouped_ms2_data)
}

#' Filter Low Complexity Spectra
#'
#' Filter peaks and their associated spectra based on the number of fragments accounting for a fraction of the signal (i.e., > 5 fragments each account for more than 1% of the signal)
#'
#' @inheritParams filter_low_information_fragments
#' @param min_signal_fragments minimum number of fragments that a spectra must have to be usable
#'
#' @return grouped_ms2_data with low-complexity spectra and their fragments removed.
filter_low_complexity_spectra <- function(grouped_ms2_data, min_signal_fragments = 5, signal_fraction_co = 0.01) {
  stopifnot("data.frame" %in% class(grouped_ms2_data))
  stopifnot(class(min_signal_fragments) == "numeric", length(min_signal_fragments) == 1, min_signal_fragments >= 0)
  stopifnot(class(signal_fraction_co) == "numeric", length(signal_fraction_co) == 1, signal_fraction_co > 0, signal_fraction_co < 1)


  complex_spectra <- filter_low_information_fragments(grouped_ms2_data, filter_method = "signal_pool", signal_fraction_co = signal_fraction_co) %>%
    dplyr::group_by(peak_num) %>%
    dplyr::filter(dplyr::n() >= min_signal_fragments) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(peak_num)

  grouped_ms2_data %>%
    dplyr::semi_join(complex_spectra, by = "peak_num")
}

plot_ms1_group <- function(mz_set, grouped_ms2_data, rt_range, label_groups = FALSE, plot_title = "") {

  # renumber mz_set
  grouped_ms2_data <- grouped_ms2_data %>%
    dplyr::left_join(grouped_ms2_data %>%
      dplyr::ungroup() %>%
      dplyr::distinct(mz_set) %>%
      dplyr::arrange(mz_set) %>%
      dplyr::mutate(mz_set_new = 1:dplyr::n()),
    by = "mz_set"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-mz_set) %>%
    dplyr::rename(mz_set = mz_set_new)

  # cluster by MS2 groups using cosine similarity
  grouped_ms2_data_hclust <- grouped_ms2_data %>%
    reshape2::acast(peak_num ~ mz_set, value.var = "ms2_ic", fill = 0) %>%
    cosine_similarity() %>%
    stats::dist() %>%
    # hclust can figure out the a similarity measure was provided rather than distance
    stats::hclust()

  # reorder MS2 groups by cosine similarity
  grouped_ms2_data$peak_num <- factor(grouped_ms2_data$peak_num, levels = grouped_ms2_data_hclust$order)

  # Prepare all of the legends so that they can be seperately created

  # Retention time and MS1 IC color scheme

  rt_colors <- tibble::tibble(rt = seq(from = rt_range[1], to = rt_range[2], length.out = 1000), color = grDevices::rainbow(1000))
  ic_colors <- tibble::tibble(log_ic = seq(from = log(100), to = log(1e8), length.out = 1000), color = grDevices::colorRampPalette(c("blue", "yellow"))(1000)) %>%
    dplyr::mutate(ic = exp(log_ic))

  peak_rt_ic <- mz_set %>%
    # select for fields allowing for an optional group
    dplyr::select_(.dots = intersect(c("peak_num", "mz_set", "precursorIc", "rt", "group"), colnames(mz_set))) %>%
    dplyr::mutate(peak_num = factor(peak_num, levels = grouped_ms2_data_hclust$order))
  peak_rt_ic$ic_color <- sapply(peak_rt_ic$precursorIc, function(an_ic) {
    ic_colors$color[which.min(abs(ic_colors$ic - an_ic))]
  })
  peak_rt_ic$rt_color <- sapply(peak_rt_ic$rt, function(an_rt) {
    rt_colors$color[which.min(abs(rt_colors$rt - an_rt))]
  })

  if (label_groups) {
    group_colors <- tibble::tibble(group = unique(mz_set$group)) %>%
      dplyr::mutate(group_color = grDevices::rainbow(dplyr::n()))
    peak_rt_ic <- peak_rt_ic %>%
      dplyr::left_join(group_colors, by = "group")
  }

  # generate MS2 signal intensity color scheme
  MS2_ic_frac_colors <- tibble::tibble(
    lms2_ic_frac = seq(from = log(1e-6), to = log(1), length.out = 100),
    color = grDevices::colorRampPalette(c("gray50", "black", "red"))(100)
  ) %>%
    dplyr::mutate(ms2_ic_frac = exp(lms2_ic_frac))

  grouped_ms2_data$color <- sapply(grouped_ms2_data$ms2_ic_frac, function(an_ic_frac) {
    MS2_ic_frac_colors$color[which.min(abs(MS2_ic_frac_colors$ms2_ic_frac - an_ic_frac))]
  })

  n_mz_set <- max(grouped_ms2_data$mz_set)

  x_labels_mz <- grouped_ms2_data %>%
    dplyr::group_by(mz_set, ms2_mz) %>%
    dplyr::summarize(ic_frac_sum = sum(ms2_ic_frac)) %>%
    dplyr::arrange(desc(ic_frac_sum)) %>%
    dplyr::ungroup() %>%
    dplyr::slice(1:10)

  # generate plots

  mz_group_theme <- theme_minimal() +
    theme(
      text = element_text(size = 20, color = "black"),
      axis.text.x = element_text(angle = 60, hjust = 1),
      axis.text.y = element_blank(),
      title = element_text(size = 28, color = "black"),
      panel.grid = element_blank()
    )

  MS2_plot <- ggplot(grouped_ms2_data, aes(y = peak_num)) +
    # MS2 parent peak ID ~ MS2 fingerprint
    geom_tile(aes(x = mz_set, fill = color)) +
    # RT of parent peaks
    geom_tile(data = peak_rt_ic, aes(x = n_mz_set * 1.01 + n_mz_set / 60, fill = rt_color), width = n_mz_set / 30) +
    # IC of parent peaks
    geom_tile(data = peak_rt_ic, aes(x = n_mz_set * 1.01 + n_mz_set / 20, fill = ic_color), width = n_mz_set / 30) +
    scale_fill_identity() +
    # scale_x_continuous("MS2: M/Z Bin", expand = c(0,0)) +
    scale_x_continuous("MS2: M/Z Bin", breaks = x_labels_mz$mz_set, labels = x_labels_mz$ms2_mz, expand = c(0, 0)) +
    scale_y_discrete("Peaks within MS1 M/Z tolerance", expand = c(0, 0)) +
    mz_group_theme +
    ggtitle(plot_title)

  if (label_groups) {
    MS2_plot <- MS2_plot +
      geom_tile(data = peak_rt_ic, aes(x = 0 - n_mz_set / 60, fill = group_color), width = n_mz_set / 30) +
      geom_text(data = data.frame(x = 0 - n_mz_set / 60, peak_num = levels(peak_rt_ic$peak_num)[ceiling(nrow(peak_rt_ic) / 2)]), aes(x = x, label = "Groups"), angle = 90, size = 5)
  }

  rt_colorbar <- make_colorbar(
    values = rt_colors$rt,
    colors = rt_colors$color,
    is_ln = FALSE,
    nbreaks = 10,
    label = "Retention\nTime"
  )

  ic_colorbar <- make_colorbar(
    values = ic_colors$ic,
    colors = ic_colors$color,
    is_ln = TRUE,
    nbreaks = 10,
    label = "Parent\nIC"
  )

  MS2_colorbar <- make_colorbar(
    values = MS2_ic_frac_colors$ms2_ic_frac,
    colors = MS2_ic_frac_colors$color,
    is_ln = TRUE,
    nbreaks = 10,
    label = "MS2\nSignal\nFraction"
  )

  gridExtra::grid.arrange(grobs = list(MS2_plot, rt_colorbar, ic_colorbar, MS2_colorbar), layout_matrix = cbind(matrix(1, nrow = 3, ncol = 7), matrix(2:4, nrow = 3)))
}

make_colorbar <- function(values, colors, is_ln, nbreaks, label) {
  colorbar_data <- tibble::tibble(x = 1, y = values, color = colors)

  three_sig <- function(x) {
    signif(x, 3)
  }

  ggplot(colorbar_data, aes(x = x, y = y, fill = colors)) +
    geom_raster() +
    scale_fill_identity() +
    scale_x_discrete(NULL, expand = c(0, 0)) +
    scale_y_continuous(label,
      trans = ifelse(is_ln, "log", "identity"),
      breaks = unname(stats::quantile(values, seq(0, 1, length.out = nbreaks))),
      labels = three_sig
    ) +
    theme(
      line = element_blank(),
      rect = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 15),
      axis.title.y = element_text(color = "black", size = 20, angle = 0, vjust = 0.5)
    )
}

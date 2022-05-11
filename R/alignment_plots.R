#' MS2 Driven Aligner Plotting
#'
#' @param rt_plot_data plot_data produced using \code{\link{ms2_driven_aligner}}.
#' @param print_plots If FALSE then return individual plots in a list; if TRUE then print a plot and silently return 1.
#' @param show_specific_samples A vector of sampleIds or if NULL a set of sampleIds will be randomely sampled. If sampleIds are named then they will be used for plots.
#'
#' @details
#' \itemize{
#'     \item{1: dendrogram comparing samples based on shared MS2 manhattan distance}
#'     \item{2: tile plot comparing samples based on shared MS2 counts}
#'     \item{3: visualize RT deviation trend for individual samples}
#'     \item{4: compare the RT of MS/MS2 matched samples for pairs of samples}
#'     \item{5: consensus ppm drift of shown samples}
#'     \item{6: mass accuracy deviations of individual samples}
#'     }
#'
#' @export
ms2_driven_aligner_plotting <- function(rt_plot_data, print_plots, show_specific_samples = NULL) {
  stopifnot(class(print_plots) == "logical", length(print_plots) == 1, print_plots %in% c(TRUE, FALSE))
  stopifnot(class(show_specific_samples) %in% c("NULL", "character", "integer", "numeric"))

  # check for mis-specified entries or choose random samples if NULL

  show_specific_samples <- validate_show_specific_samples(show_specific_samples, rt_plot_data$samples)

  # Summarize sample similarity based on co-detection of features
  codetection_plots <- ms2_driven_aligner_plotting_codetection(rt_plot_data, show_specific_samples)

  # Summarizing sample-to-sample deviations in retention time
  rt_alignment_plots <- ms2_driven_aligner_plotting_rt_alignment(rt_plot_data, show_specific_samples)

  # Summarizing sample-2-sample ppm drift
  mz_deviation_plots <- ms2_driven_aligner_plotting_mz_deviation(rt_plot_data, show_specific_samples)

  output <- c(codetection_plots, rt_alignment_plots, mz_deviation_plots)

  if (print_plots == FALSE) {
    return(output)
  }

  if (all(c("sample_codetection_tile", "rt_deviation", "rt_2_sample_comparison", "mz_ppm_bias", "mz_deviation") %in% names(output))) {
    gridExtra::grid.arrange(output$sample_codetection_tile,
      output$rt_deviation,
      output$rt_2_sample_comparison,
      output$mz_ppm_bias,
      output$mz_deviation,
      layout_matrix = rbind(
        c(1, 1, 3, 3),
        c(1, 1, 3, 3),
        c(4, 4, 2, 2),
        c(5, 5, 2, 2),
        c(5, 5, 2, 2)
      )
    )
    invisible(0)
  } else {
    stop("Some plots could not be generated so no plots will be printed; call this function with print_plots = FALSE to see individual plots")
  }
}

validate_show_specific_samples <- function(show_specific_samples, samples) {
  all_samples <- samples$sampleId

  if (!is.null(show_specific_samples)) {
    missing_specific_samples <- setdiff(show_specific_samples, all_samples)
    if (length(missing_specific_samples) == length(unique(show_specific_samples))) {
      stop("No samples present among show_specific_samples: ", paste(show_specific_samples, collapse = ", "))
    }
    if (length(missing_specific_samples) != 0) {
      warning("Samples from show_specific_samples not found: ", paste(missing_specific_samples, collapse = ", "))
    }
  } else {
    show_specific_samples <- sort(sample(all_samples, pmin(length(all_samples), 10)))
  }

  if (class(names(show_specific_samples)) == "NULL") {
    # add names since they can be provided in show_specific_samples for labelling
    names(show_specific_samples) <- as.character(show_specific_samples)
  }

  if (length(unique(names(show_specific_samples))) != length(names(show_specific_samples))) {
    stop("some names in \"show_specific_samples\" are duplicated; all names must be unique")
  }

  return(show_specific_samples)
}

ms2_driven_aligner_plotting_codetection <- function(rt_plot_data, show_specific_samples, max_shared_count_label = 200) {
  stopifnot(
    class(max_shared_count_label) %in% c("numeric", "integer"),
    length(max_shared_count_label) == 1,
    max_shared_count_label >= 0
  )

  ms2_group_detections <- rt_plot_data$sample_feature_codetection %>%
    dplyr::filter(sampleId %in% show_specific_samples) %>%
    dplyr::left_join(tibble::tibble(
      sampleId = unname(show_specific_samples),
      sampleLabel = names(show_specific_samples)
    ),
    by = "sampleId"
    )

  sample_similarity_manhattan <- ms2_group_detections %>%
    dplyr::select(-sampleId) %>%
    dplyr::mutate(present = 1) %>%
    tidyr::spread("sampleLabel", "present", fill = 0) %>%
    dplyr::select(-mz_set, -ms2_group) %>%
    as.matrix() %>%
    t() %>%
    stats::dist(method = "manhattan") %>%
    stats::hclust()

  output <- list()
  output$sample_codetection_dendrogram <- ggdendro::ggdendrogram(sample_similarity_manhattan, rotate = TRUE) +
    ggtitle("Dendrogram of all samples using manhattan distance")

  sample_codetection_pairs <- ms2_group_detections %>%
    dplyr::rename(sampleId_1 = sampleId, sampleLabel_1 = sampleLabel) %>%
    dplyr::left_join(ms2_group_detections %>%
      dplyr::rename(sampleId_2 = sampleId, sampleLabel_2 = sampleLabel),
    by = c("mz_set", "ms2_group")
    ) %>%
    dplyr::count(sampleId_1, sampleId_2, sampleLabel_1, sampleLabel_2) %>%
    dplyr::mutate(
      sampleLabel_1 = factor(sampleLabel_1, levels = names(show_specific_samples)),
      sampleLabel_2 = factor(sampleLabel_2, levels = names(show_specific_samples))
    )

  output$sample_codetection_tile <- ggplot(sample_codetection_pairs, aes(x = sampleLabel_1, y = sampleLabel_2, fill = n)) +
    geom_raster() +
    geom_text(
      data = sample_codetection_pairs %>%
        dplyr::filter(n < max_shared_count_label),
      aes(label = n), color = "black"
    ) +
    scale_fill_gradientn("Shared MS2s", colours = c("blue", "white", "red"), trans = "sqrt") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank()
    )

  return(output)
}

ms2_driven_aligner_plotting_rt_alignment <- function(rt_plot_data, show_specific_samples) {
  output <- list()

  if ("rt_overall_deviations" %in% names(rt_plot_data)) {
    selected_rt_deviations <- rt_plot_data$rt_overall_deviations %>%
      dplyr::semi_join(tibble::tibble(sampleId = show_specific_samples), by = "sampleId") %>%
      dplyr::left_join(tibble::tibble(
        sampleId = unname(show_specific_samples),
        sampleLabel = names(show_specific_samples)
      ),
      by = "sampleId"
      ) %>%
      dplyr::select(-sampleId) %>%
      dplyr::mutate(sample = factor(sampleLabel, levels = names(show_specific_samples))) %>%
      dplyr::arrange(rt)

    if (nrow(selected_rt_deviations) > 0) {
      rt_range <- range(selected_rt_deviations$rt)

      output$rt_deviation <- ggplot(selected_rt_deviations) +
        geom_point(aes(x = rt, y = rt_residual, alpha = is_retained), col = "blue", size = 0.5) +
        geom_point(aes(x = rt, y = rt_residual - rt_deviation, alpha = is_retained), col = "red", size = 0.5) +
        geom_path(aes(x = rt, y = rt_deviation)) +
        geom_hline(yintercept = 0, col = "black", alpha = 0.5) +
        scale_y_continuous("RT observed - mean across samples of deviation adjusted RT") +
        coord_cartesian(ylim = c(-1 * diff(rt_range) / 10, diff(rt_range) / 10)) +
        facet_wrap(~sample) +
        scale_alpha_manual(values = c("FALSE" = 0.1, "TRUE" = 0.3), guide = FALSE) +
        ggtitle("Sample retention time deviation", subtitle = "matched peaks (blue), deviation fit (black) and updated (red)") +
        theme_bw()
    }
  }

  if ("rt_2_sample_comparison" %in% names(rt_plot_data)) {
    sample_overlap <- rt_plot_data$rt_2_sample_comparison %>%
      dplyr::filter(
        s1 %in% show_specific_samples,
        s2 %in% show_specific_samples
      ) %>%
      # s1, s2 should be integers not ordered
      dplyr::mutate(
        s1 = as.integer(as.character(s1)),
        s2 = as.integer(as.character(s2))
      ) %>%
      dplyr::left_join(tibble::tibble(
        s1 = unname(show_specific_samples),
        s1l = names(show_specific_samples)
      ),
      by = "s1"
      ) %>%
      dplyr::left_join(tibble::tibble(
        s2 = unname(show_specific_samples),
        s2l = names(show_specific_samples)
      ),
      by = "s2"
      ) %>%
      dplyr::mutate(
        s1l = factor(s1l, levels = names(show_specific_samples)),
        s2l = factor(s2l, levels = names(show_specific_samples))
      )

    if (nrow(sample_overlap) > 0) {
      output$rt_2_sample_comparison <- ggplot(sample_overlap) +
        geom_point(aes(x = rt_1_raw, y = rt_2_raw), color = "blue", size = 0.8) +
        geom_point(aes(x = rt_1_update, y = rt_2_update), color = "red", size = 0.8) +
        geom_abline(intercept = 0, slope = 1, color = "black") +
        facet_grid(s2l ~ s1l) +
        scale_x_continuous("RT: sample x") +
        scale_y_continuous("RT: sample y") +
        ggtitle("Retention times of matched compounds between pairs of samples", subtitle = "before (blue) and after (red) RT alignment") +
        theme_bw()
    }
  }

  return(output)
}

ms2_driven_aligner_plotting_mz_deviation <- function(rt_plot_data, show_specific_samples) {
  output <- list()

  if ("mz_ppm_drift" %in% names(rt_plot_data) && nrow(rt_plot_data$mz_ppm_drift) > 0) {
    output$mz_ppm_bias <- ggplot(rt_plot_data$mz_ppm_drift, aes(x = ppm_deviation, y = ..density..)) +
      geom_histogram(binwidth = 0.05, fill = "blue") +
      theme_bw() +
      ggtitle("consensus ppm bias of all samples") +
      theme(text = element_text(size = 15))
  }

  if ("mz_ppm_deviations" %in% names(rt_plot_data)) {
    selected_mz_deviations <- rt_plot_data$mz_ppm_deviations %>%
      dplyr::semi_join(tibble::tibble(sampleId = show_specific_samples), by = "sampleId") %>%
      dplyr::left_join(tibble::tibble(
        sampleId = unname(show_specific_samples),
        sampleLabel = names(show_specific_samples)
      ),
      by = "sampleId"
      ) %>%
      dplyr::select(-sampleId) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sampleLabel = factor(sampleLabel, levels = names(show_specific_samples)))

    selected_mz_drift <- rt_plot_data$mz_ppm_drift %>%
      dplyr::semi_join(tibble::tibble(sampleId = show_specific_samples), by = "sampleId") %>%
      dplyr::left_join(tibble::tibble(
        sampleId = unname(show_specific_samples),
        sampleLabel = names(show_specific_samples)
      ),
      by = "sampleId"
      ) %>%
      dplyr::select(-sampleId) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sampleLabel = factor(sampleLabel, levels = names(show_specific_samples)))

    if (nrow(selected_mz_deviations) > 0) {
      output$mz_deviation <- ggplot(selected_mz_deviations, aes(x = ppm_deviation, y = ..density..)) +
        geom_histogram(binwidth = 0.1) +
        geom_vline(xintercept = 0, size = 1, color = "yellow") +
        geom_vline(
          data = selected_mz_drift,
          aes(xintercept = ppm_deviation), color = "chartreuse"
        ) +
        scale_x_continuous("ppm deviation") +
        facet_wrap(~sampleLabel, ncol = 2, strip.position = "right") +
        coord_cartesian(xlim = c(-5, 5)) +
        theme_bw() +
        ggtitle("ppm deviation between matched compounds of a sample and the consensus") +
        theme(text = element_text(size = 15))
    }
  }

  return(output)
}

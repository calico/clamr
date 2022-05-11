#' Label Coelutions
#'
#' @description Apply labels to peak groups based on sample-level coelutions defined by mz differences characteristic of common adducts and isotope mass differences.
#'
#' @details
#' Based on set of valid mz-differences for coelutions (a coelution object containing a set of valid coelution_mz_diffs), build clusters of mzs which are linked via adducts and isotopologues.
#' Pairs of [mz, scan] coelutions which are aggregated if either mz-value is within the supplied mass accuracy tolerance and coelutions are within small scan difference.
#' Build a sample-level map of coelution relationships and then apply this to peaks.
#'
#' @inheritParams test_mzroll_db_con_schema
#' @inheritParams test_clamr_config
#' @param minor_edge_filter_fraction peakgroup pairs which coelute less than minor_edge_filter_fraction will not be used for peakgroup aggregation"
#' @inheritParams identify_coelutions
#' @inheritParams simplify_coelution_subgraph
#'
#' @export
peakgrp_coelution_annotation <- function(mzroll_db_con, clamr_config, coelutions, minor_edge_filter_fraction = 0.01, min_edge_inconsistency_for_root = 1, root_weight_constant = 1.25) {
  debugr::dwatch(msg = "started peakgrp_coelution_annotation [clamr<peakgrp_coelution_annotation.R>::peakgrp_coelution_annotation]")

  # test mzroll_db_con and clamr_config

  test_mzroll_db_con_schema(mzroll_db_con)
  stopifnot(class(minor_edge_filter_fraction) == "numeric", minor_edge_filter_fraction > 0, minor_edge_filter_fraction < 1)
  require_tolerances(clamr_config, required_msLevels = 1L)

  # Test other inputs
  if (!("coelution" %in% class(coelutions))) {
    stop('"coelutions" must contain the "coelution" class')
  }

  if (!all(c("library_formula", "library_is_loss", "library_type") %in% colnames(coelutions$coelutions)) |
    !all(c("library_formula", "library_is_loss", "library_type") %in% colnames(coelutions$coelution_mz_diffs))) {
    stop('"library_formula", "library_is_loss", "library_type" were not all present as variables in the "coelutions"
and "coelutions_mz_diffs" table of the coelution object: add these fields with the "apply_coelution_labels" function')
  }

  # Determine whether the same samples are provided in the clamr_db as the coelutions

  # load the sample table mzroll_db
  all_mzroll_samples <- dplyr::tbl(mzroll_db_con, "samples") %>%
    dplyr::collect() %>%
    dplyr::mutate(name_no_extn = strip_ms_file_extension(name))
  # summarize the unique samples in coelutions
  all_coelution_samples <- unique(coelutions$coelutions$sample)

  debugr::dwatch(msg = "retrieved samples information from mzrollDB. [clamr<peakgrp_coelution_annotation.R>::peakgrp_coelution_annotation]")

  mismatched_samples <- setdiff(union(all_mzroll_samples$name_no_extn, all_coelution_samples), intersect(all_mzroll_samples$name_no_extn, all_coelution_samples))
  if (length(mismatched_samples) != 0) {
    warning(length(mismatched_samples), " samples are not matched between the coelutions summary and the clamr_db")
    warning(length(setdiff(all_mzroll_samples$name_no_extn, all_coelution_samples)), " are present in clamr_db but not coelutions: ", paste(setdiff(all_mzroll_samples$name_no_extn, all_coelution_samples), collapse = ", "))
    warning(length(setdiff(all_coelution_samples, all_mzroll_samples$name_no_extn)), " are present in coelutions but not clamr_db: ", paste(setdiff(all_coelution_samples, all_mzroll_samples$name_no_extn), collapse = ", "))
  }

  # Build a sample-level map of coelution relationships and then apply this to peaks

  clamr_db_peaks <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT peakId, groupId, sampleId, peakMz, mzmin, mzmax, scan, minscan, maxscan FROM peaks")) %>%
    dplyr::collect()

  debugr::dwatch(msg = "retrieved peaks information from mzrollDB. [clamr<peakgrp_coelution_annotation.R>::peakgrp_coelution_annotation]")

  # match peaks in clamr_db with peaks in coelutions
  # doing this we can pass scan-level coelutions between peaks to coelutions among peaks in a peak groups

  sample_peak_coelutions <- all_mzroll_samples %>%
    # add peaks nested within samples
    dplyr::inner_join(
      tidyr::nest(clamr_db_peaks, peaks = -sampleId),
      by = "sampleId"
    ) %>%
    # add coelutions nested within samples
    dplyr::inner_join(
      tidyr::nest(coelutions$coelutions, coelutions = -sample),
      by = c("name_no_extn" = "sample")
    ) %>%
    # for each sample, determine coelution peaks [mz, scans] and then match these to clamr db peaks
    dplyr::mutate(matches = furrr::future_map2(peaks, coelutions, combine_sample_coelutions_w_peaks, clamr_config$MS1tol)) %>%
    dplyr::select(-peaks, -coelutions) %>%
    tidyr::unnest(matches)

  debugr::dwatch(msg = "matched peaks in db with peaks in coelutions. [clamr<peakgrp_coelution_annotation.R>::peakgrp_coelution_annotation]")

  # Handle case of no coelutions detected
  if (nrow(sample_peak_coelutions) == 0) {
    debugr::dwatch(msg = "no coelutions were detected. [clamr<peakgrp_coelution_annotation.R>::peakgrp_coelution_annotation]")
    return(NULL)
  }

  group_level_coelution <- sample_peak_coelutions %>%
    # collapse 1-many clamr peak - coelution peak into unique entries (these may exist if different parts of a peak pair are separately extracted but not fully combined)
    dplyr::distinct(sampleId, peak_1_peakId, peak_1_groupId, peak_2_peakId, peak_2_groupId, library_formula, library_is_loss, library_type) %>%
    dplyr::count(peak_1_groupId, peak_2_groupId, library_formula, library_is_loss, library_type)

  # find totally separable groups
  minor_edge_co <- floor(length(all_coelution_samples) * minor_edge_filter_fraction)

  coelution_preliminary_groups <- group_level_coelution %>%
    dplyr::filter(n > minor_edge_co) %>%
    igraph::graph_from_data_frame(directed = FALSE) %>%
    igraph::clusters() %>%
    igraph::membership() %>%
    {
      tibble::tibble(groupId = as.integer(names(.)), prelim_metaGroupId = as.integer(unname(.)))
    }

  debugr::dwatch(msg = "found totally separable coelution groups. [clamr<peakgrp_coelution_annotation.R>::peakgrp_coelution_annotation]")

  coelution_subgraphs <- group_level_coelution %>%
    # match peakgroups to coelution subsets
    dplyr::inner_join(coelution_preliminary_groups, by = c("peak_1_groupId" = "groupId")) %>%
    # remove edges not contained within subset
    dplyr::semi_join(coelution_preliminary_groups, by = c("peak_2_groupId" = "groupId", "prelim_metaGroupId")) %>%
    tidyr::nest(data = -prelim_metaGroupId) %>%
    # not actually needed once peakdetector is debugged
    dplyr::rowwise() %>%
    dplyr::mutate(nrow = nrow(data))

  max_group_size <- max(coelution_subgraphs$nrow)

  if (any(coelution_subgraphs$nrow) > 5000) {
    warning("max group size is ", max_group_size, "; a group size over 5,000 may take considerable prime to process; increase")
  }

  # Among a set of coeluting peaks determine which are the monoisotopic/non-adduct peaks and define metaGroupds based on these parentGroupIds

  metaGroups <- coelution_subgraphs %>%
    {
      furrr::future_map(.$data, simplify_coelution_subgraph, min_edge_inconsistency_for_root = min_edge_inconsistency_for_root, root_weight_constant = root_weight_constant)
    } %>%
    dplyr::bind_rows()

  debugr::dwatch(msg = "found coelution subgroups metagroups. [clamr<peakgrp_coelution_annotation.R>::peakgrp_coelution_annotation]")

  # add disconnected peakGroups

  peakGroup_update <- metaGroups %>%
    dplyr::bind_rows(
      clamr_db_peaks %>%
        dplyr::distinct(groupId) %>%
        dplyr::anti_join(metaGroups, by = "groupId") %>%
        dplyr::mutate(parentGroupId = groupId, adductName = "")
    )

  peakGroup_update <- peakGroup_update %>%
    dplyr::left_join(
      peakGroup_update %>%
        dplyr::ungroup() %>%
        dplyr::distinct(parentGroupId) %>%
        dplyr::mutate(metaGroupId = 1:dplyr::n()),
      by = "parentGroupId"
    ) %>%
    # parents are assigned a parentGroupId = 0 by convention
    dplyr::mutate(parentGroupId = ifelse(parentGroupId == groupId, 0L, parentGroupId)) %>%
    dplyr::ungroup()

  debugr::dwatch(msg = "finished updating peak groups with coelution annotation information.\n[clamr<peakgrp_coelution_annotation.R>::peakgrp_coelution_annotation]")

  return(peakGroup_update)
}

#' Simplify Coelution Clique
#'
#' To encourage paths which link to a common root, among the shortest path, the path to the best root is chosen by taking the smallest path weight + root weight.
#'
#' @param one_coelution_subgraph the edgelist for a single prelim_metaGroup
#' @param min_edge_inconsistency_for_root When determining the network's roots, ignore edges represented by less than \code{min_edge_inconsistency_for_root} samples.
#' @param root_weight_constant W(R) = C / sum(n of daughter edges + 1)
#'
#' @return a tibble containing groupId summaries
simplify_coelution_subgraph <- function(one_coelution_subgraph, min_edge_inconsistency_for_root = 1, root_weight_constant = 1.25) {

  # separate a clique into shortest paths from each leaf to a root
  # convert to character since this is how igraph will store vertex names
  all_vertices <- as.character(unique(c(one_coelution_subgraph$peak_1_groupId, one_coelution_subgraph$peak_2_groupId)))

  strongly_connected_subset <- one_coelution_subgraph %>%
    dplyr::filter(n > min_edge_inconsistency_for_root)

  # define the set of vertices that could be roots of the network (no in edge or weak in edges (based on min_edge_consistency_for_rtoot))
  # a possible root does not have any parents (among the strongly connected subset)

  possible_roots <- if (nrow(strongly_connected_subset) == 0) {
    all_vertices
  } else {
    strongly_connected_subset %>%
      igraph::graph_from_data_frame(directed = TRUE) %$%
      sapply(igraph::V(.), function(a_vertex) {
        igraph::subcomponent(graph = ., a_vertex, "in")
      }) %>%
      purrr::map_dbl(length) %>%
      {
        names(.)[. == min(.)]
      } %>%
      # add nodes which were not in the strongly_connected_subset
      {
        c(., setdiff(all_vertices, as.character(unique(c(strongly_connected_subset$peak_1_groupId, strongly_connected_subset$peak_2_groupId)))))
      }
  }

  # calculate the weighted shortest path between all root and all vertices
  # A shortest path is calculated by minimizing the W(P) = sum(W(E)). W(E) = 1 / n samples w/ edge
  # To choose among paths involving different roots, a root weight, W(R), is added to W(P). W(R) = C / sum(n of daughter edges).

  # setup network, add edge weights
  full_directed_network <- one_coelution_subgraph %>%
    igraph::graph_from_data_frame(directed = TRUE)
  # weight each edge by 1 / n samples w/ edge
  igraph::E(full_directed_network)$weight <- 1 / igraph::E(full_directed_network)$n

  # find the path with the lowest summed edge weight between each root and all reachable vertices

  all_best_paths <- lapply(possible_roots, function(a_root) {
    tibble::tibble(
      root = a_root,
      vertex = all_vertices,
      epath = suppressWarnings(igraph::shortest_paths(full_directed_network,
        from = a_root,
        to = all_vertices,
        mode = "out",
        output = "epath"
      ))$epath
    )
  }) %>%
    dplyr::bind_rows() %>%
    # remove infeasible paths
    dplyr::mutate(nsteps = purrr::map_int(.$epath, length)) %>%
    dplyr::filter(nsteps != 0) %>%
    # convert from igraph epaths to tibble
    dplyr::mutate(path = purrr::map(epath, expand_epath)) %>%
    tidyr::unnest(path) %>%
    tidyr::separate(edge, into = c("peak_1_groupId", "peak_2_groupId")) %>%
    dplyr::mutate_at(.vars = dplyr::vars(peak_1_groupId, peak_2_groupId), .funs = dplyr::funs(as.integer)) %>%
    dplyr::left_join(one_coelution_subgraph, by = c("peak_1_groupId", "peak_2_groupId"))

  root_weights <- tibble::tibble(root = possible_roots) %>%
    dplyr::left_join(one_coelution_subgraph %>%
      dplyr::filter(peak_1_groupId %in% possible_roots) %>%
      dplyr::group_by(peak_1_groupId) %>%
      dplyr::summarize(n_daughters = sum(n)) %>%
      dplyr::mutate(peak_1_groupId = as.character(peak_1_groupId)),
    by = c("root" = "peak_1_groupId")
    ) %>%
    dplyr::mutate(
      n_daughters = ifelse(is.na(n_daughters), 0, n_daughters),
      root_weight = root_weight_constant / (n_daughters + 1)
    )

  # find the path with the lowest path weight + root weight for each vertex

  optimal_paths <- all_best_paths %>%
    dplyr::group_by(root, vertex) %>%
    dplyr::summarize(path_weight = sum(1 / n)) %>%
    # add each possible root as its own root: "I am Root"
    dplyr::bind_rows(tibble::tibble(root = possible_roots, vertex = possible_roots, path_weight = 0)) %>%
    # add root weight
    dplyr::left_join(root_weights, by = "root") %>%
    dplyr::mutate(total_path_weight = path_weight + root_weight) %>%
    # take the best path for each vertex
    dplyr::group_by(vertex) %>%
    dplyr::arrange(path_weight + root_weight) %>%
    dplyr::slice(1)

  # write out the parent ion (root) for each peak group and the coelution path

  parentAssignments <- optimal_paths %>%
    dplyr::select(root, vertex) %>%
    dplyr::left_join(all_best_paths, by = c("root", "vertex")) %>%
    dplyr::mutate(
      groupId = as.integer(vertex),
      parentGroupId = as.integer(root)
    ) %>%
    dplyr::group_by(groupId, parentGroupId) %>%
    dplyr::summarize(adductName = paste(library_formula, collapse = ", ")) %>%
    dplyr::mutate(adductName = ifelse(adductName == "NA", "", adductName))

  # check that all vertices were assigned a parent (even if they are their own parent)
  stopifnot(length(all_vertices) == nrow(parentAssignments))
  stopifnot(length(setdiff(as.integer(all_vertices), parentAssignments$groupId)) == 0)

  parentAssignments
}

expand_epath <- function(epath) {
  attr(epath, "vnames") %>%
    {
      tibble::tibble(step = seq_along(.), edge = .)
    }
}

plot_coelution_subgraph <- function(one_coelution_subgraph, group_mzs) {
  group_mz_subset <- group_mzs %>%
    dplyr::filter(groupId %in% c(one_coelution_subgraph$peak_1_groupId, one_coelution_subgraph$peak_2_groupId))

  igraph_viz <- one_coelution_subgraph %>%
    igraph::graph_from_data_frame(
      directed = TRUE,
      vertices = group_mz_subset
    )
  igraph::V(igraph_viz)$label <- igraph::V(igraph_viz)$mz
  igraph::V(igraph_viz)$size <- 10

  igraph::E(igraph_viz)$label <- paste0(igraph::E(igraph_viz)$library_formula, " (n = ", igraph::E(igraph_viz)$n, ")")
  igraph::E(igraph_viz)$width <- igraph::E(igraph_viz)$n * 0.2
  igraph::E(igraph_viz)$color <- dplyr::case_when(
    one_coelution_subgraph$library_type == "isotopologue" ~ "GREEN",
    one_coelution_subgraph$library_type == "adduct" ~ "BLUE",
    one_coelution_subgraph$library_type == "distinct compound" ~ "gray50"
  )

  if (all(c("adductName", "metaGroupId") %in% colnames(group_mzs))) {
    # bonus aesthetics

    color_palette <- group_mz_subset %>%
      dplyr::distinct(metaGroupId) %>%
      dplyr::mutate(color = grDevices::rainbow(dplyr::n()))

    group_mz_subset <- group_mz_subset %>%
      dplyr::left_join(color_palette, by = "metaGroupId")

    igraph::V(igraph_viz)$color <- group_mz_subset$color
    igraph::V(igraph_viz)$shape <- ifelse(group_mz_subset$adductName == "", "square", "circle")
  }

  igraph::plot.igraph(igraph_viz, layout = igraph::layout_with_kk, edge.arrow.size = 0.25, vertex.color = igraph::V(igraph_viz)$color)
}


#' Generate Coelution Peaks
#'
#' @description Using a set of coelutions between pairs of mzs aggregate these mzs acrosss all of the coelutions that they are involved in to
#' create a summary of peak locations [mz, scans]
#'
#' @param sample_coelutions A list generated by \code{find_common_coelutions} and possibly reduced using \code{apply_coelution_labels}.
#' @param MS1tol from \code{\link{build_clamr_config}} or \code{\link{format_mass_accuracy_input}}
#'
#' @return A list containing:
#' \itemize{
#'   \item{binned_sample_coelutions: coelutions between unique coeluting peaks}
#'   \item{unique_coeluting_peaks: unique mz-scan groups for a given sample}
#'   }
#'
#' @export
generate_coelution_peaks <- function(sample_coelutions, MS1tol) {
  coeluting_peaks <- sample_coelutions %>%
    dplyr::select(peak_1, peak_2, start_scan, final_scan) %>%
    dplyr::mutate(corr_length = final_scan - start_scan) %>%
    dplyr::mutate(row_number = 1:dplyr::n()) %>%
    tidyr::gather(mz_number, mz, peak_1, peak_2) %>%
    dplyr::arrange(mz) %>%
    # find distincts m/z up to instrument tolerance
    dplyr::mutate(mz_node_mz_set = find_mz_jumps(.$mz, MS1tol$tol, MS1tol$absolute_or_relative, collapse = FALSE)) %>%
    # group mzs into distinct peaks based on scan position
    dplyr::group_by(mz_node_mz_set) %>%
    dplyr::arrange(final_scan) %>%
    # test for whether for a given mzs scans are nearby:
    # scan[i] + tolerance[i] + tolerance[i + 1] > scan[i + 1]
    # tolerance is set as corr_length (# of scans of coelution) / 4
    # so, EIC pairs will be aggregated if they overlap for half of the EIC range
    dplyr::mutate(mz_node_scan_set = find_scan_jumps(final_scan, corr_length / 3)) %>%
    tidyr::unite(mz_scan_group, mz_node_mz_set, mz_node_scan_set)

  # define the coelution edges based on the above groups

  binned_sample_coelutions <- sample_coelutions %>%
    dplyr::bind_cols(coeluting_peaks %>%
      dplyr::select(row_number, mz_number, mz_scan_group) %>%
      tidyr::spread(mz_number, mz_scan_group) %>%
      dplyr::rename(peak_1_group = peak_1, peak_2_group = peak_2))

  # obtain a summary of a single peak with a characteristic [mz,scan]

  unique_coeluting_peaks <- coeluting_peaks %>%
    dplyr::group_by(mz_scan_group) %>%
    dplyr::summarize(
      minscan = min(final_scan - corr_length),
      maxscan = max(final_scan),
      peakMz = mean(mz),
      mzmin = min(mz),
      mzmax = max(mz),
      n = dplyr::n()
    )

  list(
    binned_sample_coelutions = binned_sample_coelutions,
    unique_coeluting_peaks = unique_coeluting_peaks
  )
}


#' Combine Sample Coelutions with Peaks
#'
#' For a single sample, from coelution mz pairs (essentially edges) infer the peaks [mz, scans] then join these coelution peaks with peaks in the clamr_db.
#'
#' @param sample_peaks peaks from a single sample from .mzrollDB
#' @param sample_coelutions coelutions from a single sample from .mzrollDB
#' @inheritParams generate_coelution_peaks
#'
#' @details
#' The edge-list defines the connection between pairs of mzs over a scan interval.
#' In order to build a graph from these connections, we need to be able to define repeated observation of the same node in terms of a common mz and scan.
combine_sample_coelutions_w_peaks <- function(sample_peaks, sample_coelutions, MS1tol) {

  # identify unique [mz,scan] groups in coelutions

  coelution_peaks <- generate_coelution_peaks(sample_coelutions, MS1tol)

  # search peaks from the clamr_db are associated with coelution peaks

  peak_to_peak_matches <- mz_peak_join_fxn(sample_peaks, coelution_peaks$unique_coeluting_peaks) %>%
    # make sure that coelution peaks are only assigned to a single clamr_db peak
    # clamr_db peaks can have multiple coelution peaks assigned
    dplyr::group_by(mz_scan_group) %>%
    dplyr::arrange(abs(peakMz.x - peakMz.y), abs((minscan.x + maxscan.x) / 2 - (minscan.y + maxscan.y) / 2)) %>%
    dplyr::slice(1) %>%
    dplyr::select(peakId, groupId, mz_scan_group)

  coelution_peaks$binned_sample_coelutions %>%
    dplyr::inner_join(peak_to_peak_matches %>%
      dplyr::rename(peak_1_group = mz_scan_group, peak_1_peakId = peakId, peak_1_groupId = groupId),
    by = "peak_1_group"
    ) %>%
    dplyr::inner_join(peak_to_peak_matches %>%
      dplyr::rename(peak_2_group = mz_scan_group, peak_2_peakId = peakId, peak_2_groupId = groupId),
    by = "peak_2_group"
    )
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peakgrp_coelution_annotation.R
\name{simplify_coelution_subgraph}
\alias{simplify_coelution_subgraph}
\title{Simplify Coelution Clique}
\usage{
simplify_coelution_subgraph(
  one_coelution_subgraph,
  min_edge_inconsistency_for_root = 1,
  root_weight_constant = 1.25
)
}
\arguments{
\item{one_coelution_subgraph}{the edgelist for a single prelim_metaGroup}

\item{min_edge_inconsistency_for_root}{When determining the network's roots, ignore edges represented by less than \code{min_edge_inconsistency_for_root} samples.}

\item{root_weight_constant}{W(R) = C / sum(n of daughter edges + 1)}
}
\value{
a tibble containing groupId summaries
}
\description{
To encourage paths which link to a common root, among the shortest path, the path to the best root is chosen by taking the smallest path weight + root weight.
}

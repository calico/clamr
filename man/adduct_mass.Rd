% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemistry_utils.R
\name{adduct_mass}
\alias{adduct_mass}
\title{Adduct Mass}
\usage{
adduct_mass(adducts)
}
\arguments{
\item{adducts}{a character vector of adduct formulas}
}
\value{
a tibble containing one row per unique adduct
}
\description{
From a set of standard adduct formulas, e.g., [M+H]-, calculate the adducts mass and charge.
}
\examples{
adduct_mass(adducts = c("[M+H]+", "[M-NH4]2-", "[M]+"))

}

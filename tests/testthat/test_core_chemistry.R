context("Test Core Chemistry Utilities")

test_that("Formulas are properly read", {

  # formulas containing gains and losses
  expect_equal(split_formula("C6H12O6")$gains, "C6H12O6")
  expect_equal(split_formula("C6H12O6 + Na")$gains, "C6H12O6Na")
  expect_equal(split_formula("C6H12O6 - H"), list(gains = "C6H12O6", losses = "H"))
  expect_equal(split_formula("-H2")$losses, "H2")
  expect_equal(split_formula("-ZnS5 + PtO4"), list(gains = "PtO4", losses = "ZnS5"))

  expect_equal(split_formula("13C C(5) H(12) O(6) - H(5)", collapse = " "), list(gains = "13C C(5) H(12) O(6)", losses = "H(5)"))

  # standard formulas
  expect_equal(parse_standard_formulas("C6H12O6")$n, c(6, 12, 6))
  # invalid formulas
  expect_warning(parse_standard_formulas("C6 H12 O6"))
  expect_warning(parse_standard_formulas("C6H12O6-H"))

  # simple formulas
  expect_equal(parse_simple_formulas("6C")$n, 6)
  expect_equal(parse_simple_formulas("H")$n, 1)
  expect_equal(parse_simple_formulas("5Na")$element, "Na")

  # isotopic formulas
  isotopic_formulas <- parse_isotopic_chemical_formula("13C C(5) H(12) O(6) - H(2)")
  expect_equal(isotopic_formulas$count, c(1, 5, 10, 6))

  # adduct formulas
  expect_equal(parse_adduct_changes("-2H")$n, -2)
  expect_equal(parse_adduct_changes("+NH4-H")$n, c(3, 1))
  
  adducts <- c("[M-2H]2-", "[M-2H+3Na]+", "[M-H]-", "[M+H]+", "[M+HCO2]-", "[M+NH4]+", "[M]+")
  expect_equal(nrow(parse_adducts(adducts)), 17)
  
  adduct_masses <- adduct_mass(adducts)
  expect_equal(nrow(adduct_masses), length(adducts))
  adduct_masses <- adduct_masses %>%
    # order so that attributes will be in the expected order
    dplyr::mutate(adduct = factor(adduct, levels = adducts)) %>%
    dplyr::arrange(adduct)
  
  expect_equal(
    adduct_masses$adduct_mass,
    c(-2.014553, 66.953109, -1.007276, 1.007276, 44.998203, 18.033826, -0.00054858),
    tolerance = 0.001
    )
  
  expect_equal(adduct_masses$adduct_charge, c(-2, 1, -1, 1, -1, 1, 1))
  
  # parse charges
  expect_equal(parse_charge_formulas(c("2-", "3+", "+"))$n, c(2, -3, -1))

  # adduct masses
  expect_equal(
    adduct_mass(adducts = "[M+H]+"),
    tibble::tibble(adduct = "[M+H]+", adduct_mass = 1.00727, adduct_charge = 1L),
    tolerance = 1e-5
    )
  expect_equal(
    adduct_mass(adducts = "[M-NH4]2-"),
    tibble::tibble(adduct = "[M-NH4]2-", adduct_mass = -18.03327, adduct_charge = -2L),
    tolerance = 1e-5
    )

})


test_that("Molecular formulas generate accurate exact masses", {

  formulas = tibble::tibble(formula = c("C6H12O6", "C5H10N2O3", "C10H17N3O6S"),
                            true_mass = c(180.063388, 146.069142, 307.083806))

  monoisotopic_masses <- formula_monoisotopic_mass(formulas$formula)

  comparison_masses <- formulas %>%
    dplyr::left_join(monoisotopic_masses, by = "formula")

  expect_equal(comparison_masses$true_mass, comparison_masses$exact_mass)

})



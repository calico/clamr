context("Test searching for standard-query matches")

test_that("searching with MZ works (one variable)", {

  query_dataset = tibble::tibble(query_id = 1:3, amu = 1:3)
  standard_dataset = clamr::isotope_summaries %>% dplyr::select(z, label, amu)
  variable_tolerances = tibble::tibble(variable = "amu", tolerance = 0.01, relative_or_absolute = "absolute")
  mz_match <- join_matching_standards(query_dataset, standard_dataset, variable_tolerances, distance_method = "manhattan", threshold = 1)

  expect_equal(nrow(mz_match), 1)
  expect_true(mz_match$query_id == 1)
  expect_true(mz_match$label == "H")
  expect_equal(round(mz_match$match_distance, 3), 0.783)

  query_dataset = tibble::tibble(query_id = 1:3, amu = c(1, 1.01, 1.02))
  mz_match <- join_matching_standards(query_dataset, standard_dataset, variable_tolerances, distance_method = "manhattan", threshold = 1)
  expect_equal(nrow(mz_match), 2)

})

test_that("Test searching with MZ + RT works (2+ variables)", {

  query_dataset = tibble::tibble(query_id = 1:3, amu = 1:3, rt = c(10, 15, 3))
  standard_dataset = clamr::isotope_summaries %>% dplyr::select(z, label, amu) %>% dplyr::mutate(rt = 1:dplyr::n())
  variable_tolerances = tibble::tribble(~variable, ~tolerance, ~relative_or_absolute,
                                        "amu", 0.01, "absolute",
                                        "rt", 5, "absolute")

  mz_match <- join_matching_standards(query_dataset, standard_dataset, variable_tolerances, distance_method = "euclidean", threshold = 1)

  expect_null(mz_match)

  mz_match <- join_matching_standards(query_dataset, standard_dataset, variable_tolerances, distance_method = "euclidean", threshold = 4)

  expect_equal(nrow(mz_match), 3)
  expect_true(all(mz_match$query_id == 1:3))
  expect_true(all(mz_match$label == c("H", "D", "3He")))
  expect_equal(round(mz_match$match_distance, 3), c(1.781, 2.784, 1.615))
})

test_that("Test searching with a relative error: PPM", {

  query_dataset = tibble::tibble(query_id = 1:4, mz = c(1, 10, 100, 1000), rt = 1:4)
  standard_dataset = tibble::tibble(library_id = 1:4, mz = c(1.001, 10.001, 100.001, 1000.001), rt = 1:4 + 0.5)
  variable_tolerances = tibble::tribble(~variable, ~tolerance, ~relative_or_absolute,
                                            "mz", 10e-6, "relative",
                                            "rt", 1, "absolute")
  mz_match <- join_matching_standards(query_dataset, standard_dataset, variable_tolerances, distance_method = "manhattan", threshold = 10)

  expect_equal(nrow(mz_match), 2)
  expect_true(all(mz_match$query_id == 3:4))
  expect_true(all(mz_match$library_id == 3:4))
  expect_equal(mz_match$match_distance, c(1.5, 0.6))
})


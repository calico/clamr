context("find_mz_jumps groups masses appropriately")

test_that("Behavior when n_max_gap_inconsistencies = 0", {

  mock_data <- c(100, 101, 101.1, 102, 200, 201)

  expect_equal(
    find_mz_jumps(mock_data, mass_accuracy_tolerance = 1, absolute_or_relative = "absolute"),
    tibble::tribble(~mz, ~mz_set,
                    100.0, 1L,
                    101.0, 2L,
                    101.1, 2L,
                    102.0, 2L,
                    200.0, 3L,
                    201.0, 4L))

  expect_equal(
    find_mz_jumps(mock_data, mass_accuracy_tolerance = 1, absolute_or_relative = "absolute", collapse = FALSE),
    c(1, 2, 2, 2, 3, 4))

  expect_equal(
    find_mz_jumps(mock_data, mass_accuracy_tolerance = 0.001, absolute_or_relative = "relative"),
    tibble::tribble(~mz, ~mz_set,
                    100.0, 1L,
                    101.0, 2L,
                    101.1, 2L,
                    102.0, 3L,
                    200.0, 4L,
                    201.0, 5L))
})

test_that("Behavior when n_max_gap_inconsistencies != 0", {

  mock_data <- c(100, 100, 100.1, 101, 101, 200, 200, 200.1, 201, 201)

  expect_equal(
    find_mz_jumps(mock_data, mass_accuracy_tolerance = 1, absolute_or_relative = "absolute", n_max_gap_inconsistencies = 1L),
    tibble::tribble(~mz, ~mz_set,
                    100.0, 1L,
                    100.1, 1L,
                    101.0, 2L,
                    200.0, 3L,
                    200.1, 3L,
                    201.0, 4L))

  expect_equal(
    find_mz_jumps(mock_data, mass_accuracy_tolerance = 0.009, absolute_or_relative = "relative", n_max_gap_inconsistencies = 1L),
    tibble::tribble(~mz, ~mz_set,
                    100.0, 1L,
                    100.1, 1L,
                    101.0, 2L,
                    200.0, 3L,
                    200.1, 3L,
                    201.0, 3L))
})

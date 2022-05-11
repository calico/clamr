context("test-Rcpp_frag_matching_algorithm")

test_that("frag matching is mz-greedy", {

  mock_experimental <- tibble::tribble(~mz, ~intensity,
                                       100.1, 1000L,
                                       200.1, 1000L)

  mock_library <- tibble::tribble(~mz, ~intensity,
                                  100.0, 1000L,
                                  101.0, 1000L,
                                  200.0, 1000L,
                                  201.0, 1000L)
  mzTol <- 5;

  matches <- mzkitcpp::msms_comparer_matches(mock_library$mz, mock_experimental$mz, mzTol, debug=TRUE)

  expect_equal(matches[1], 1) # match 100.0 <--> 100.1
  expect_equal(matches[2], -1) # no match to 101.0
  expect_equal(matches[3], 2) # match 200.1 <--> 200.0
  expect_equal(matches[4], -1) # no match to 201.0

  mock_library <- tibble::tribble(~mz, ~intensity,
                                  100.0, 1000L,
                                  101.0, 1000L,
                                  102.0, 1000L,
                                  103.0, 1000L)

  mock_experimental <- tibble::tribble(~mz, ~intensity,
                                       100.1, 1000L,
                                       101.1, 1000L,
                                       101.5, 1000L,
                                       101.9, 1000L,
                                       103.1, 1000L)

  matches <- mzkitcpp::msms_comparer_matches(mock_library$mz, mock_experimental$mz, mzTol, debug=TRUE)

  expect_equal(matches[1], 1) # match 100.0 <--> 100.1
  expect_equal(matches[2], 2) # match 101.0 <--> 101.1
  expect_equal(matches[3], 4) # match 102.0 <--> 101.9
  expect_equal(matches[4], 5) # match 103.0 <--> 103.1
})

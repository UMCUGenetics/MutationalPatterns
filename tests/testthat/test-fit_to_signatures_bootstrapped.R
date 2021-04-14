context("test-fit_to_signatures_bootstrapped")

# Get mut_mat
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

# Get signatures
signatures <- get_known_signatures()

test_that("Output has correct class", {
  output <- fit_to_signatures_bootstrapped(mut_mat, signatures, n_boots = 2, max_delta = 0.05)
  expect_true(inherits(output, "matrix"))

  output_ori <- fit_to_signatures_bootstrapped(mut_mat, signatures, n_boots = 2, max_delta = 0.05, method = "regular")
  expect_true(inherits(output_ori, "matrix"))

  output_ori_10 <- fit_to_signatures_bootstrapped(mut_mat, signatures, n_boots = 2, max_delta = 0.05, method = "regular_10+")
  expect_true(inherits(output_ori_10, "matrix"))
  
  output_best <- fit_to_signatures_bootstrapped(mut_mat, signatures[,1:5], n_boots = 2, max_delta = 0.05, method = "strict_best_subset")
  expect_true(inherits(output_best, "matrix"))
})

expected <- readRDS(system.file("states/bootstrapped_snv_refit.rds",
  package = "MutationalPatterns"
))

test_that("Output is equal to expected", {
  set.seed(42)
  output <- fit_to_signatures_bootstrapped(mut_mat, signatures, n_boots = 2, max_delta = 0.05)
  expect_equal(output, expected)
})

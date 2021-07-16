context("test-fit_to_signatures_strict")

# Get mut_mat
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

# Get signatures
signatures <- get_known_signatures()

output <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.05)
output_best <- fit_to_signatures_strict(mut_mat, signatures[,1:5], max_delta = 0.004, method = "best_subset")
output_single_sig = fit_to_signatures_strict(mut_mat, signatures[,1, drop = F], max_delta = 0.05)

expected <- readRDS(system.file("states/strict_snv_refit.rds",
  package = "MutationalPatterns"
))
expected_best <- readRDS(system.file("states/strict_best_snv_refit.rds",
                                package = "MutationalPatterns"
))

test_that("Output has correct class", {
  expect_true(inherits(output, "list"))
  expect_true(inherits(output$fit_res, "list"))
  expect_true(inherits(output$fit_res$contribution, "matrix"))
  expect_true(inherits(output$fit_res$reconstructed, "matrix"))
  expect_true(inherits(output$sim_decay_fig, "list"))
  expect_true(inherits(output$sim_decay_fig[[1]], "gg"))
  expect_true(inherits(output_best, "list"))
  expect_true(inherits(output_best$fit_res, "list"))
  expect_true(inherits(output_best$fit_res$contribution, "matrix"))
  expect_true(inherits(output_best$fit_res$reconstructed, "matrix"))
  expect_true(inherits(output_best$sim_decay_fig, "list"))
  expect_true(inherits(output_best$sim_decay_fig[[1]], "gg"))
  expect_true(inherits(output_single_sig, "list"))
  expect_true(inherits(output_single_sig$fit_res, "list"))
  expect_true(inherits(output_single_sig$fit_res$contribution, "matrix"))
  expect_true(inherits(output_single_sig$fit_res$reconstructed, "matrix"))
  expect_true(inherits(output_single_sig$sim_decay_fig, "list"))
  expect_true(inherits(output_single_sig$sim_decay_fig[[1]], "gg"))
})

test_that("Output is equal to expected", {
  expect_equal(output$fit_res, expected)
  expect_equal(output_best$fit_res, expected_best)
})

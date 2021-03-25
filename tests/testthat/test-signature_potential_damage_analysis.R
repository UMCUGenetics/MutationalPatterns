context("test-signature_potential_damage_analysis")

# Get the signatures
signatures <- get_known_signatures()

# Get the contexts
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

contexts <- rownames(mut_mat)[1:6]

# Get context mismatches
context_mismatches <- readRDS(system.file("states/context_mismatches.rds",
  package = "MutationalPatterns"
))

# Run function
output <- signature_potential_damage_analysis(signatures, contexts, context_mismatches)

test_that("Output has correct class", {
  expect_true(inherits(output, "tbl_df"))
})

test_that("Output has correct size", {
  expect_equal(dim(output), c(240, 7))
})

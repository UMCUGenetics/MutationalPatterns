context("test-cos_sim_matrix")

# Read signatures
signatures <- get_known_signatures()


# Read mut_matrix
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))


# Calculate the cosine similarity between each COSMIC signature and each 96 mutational profile
output <- cos_sim_matrix(mut_mat, signatures)

test_that("Output has correct class and data type", {
  expect_true(inherits(output, c("matrix")))
  expect_equal(typeof(output), "double")
})

test_that("Output has expected size", {
  expect_equal(dim(output), c(9, 60))
})

context("test-cos_sim_matrix")

# Read signatures
signatures <- get_known_signatures()


# Read mut_matrix
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))


# Calculate the cosine similarity between each COSMIC signature and each 96 mutational profile
output <- cos_sim_matrix(mut_mat, signatures)

# Perform tests.
test_that("Output has correct class and data type", {
  expect_true(inherits(output, c("matrix")))
  expect_equal(typeof(output), "double")
})

test_that("Output has expected size", {
  expect_equal(dim(output), c(9, 60))
})

mut_mat_df = as.data.frame(mut_mat)
mut_mat_chr = mut_mat_df
mut_mat_chr[,1] <- as.character(mut_mat_chr[,1, drop = TRUE])
test_that("Non-numeric inputs give an error", {
    expect_error(cos_sim_matrix(mut_mat_chr, signatures))
})

test_that("Tibble inputs are converted into data.frames.", {
    output2 = cos_sim_matrix(tibble::as_tibble(mut_mat_df), signatures)
    expect_equal(output, output2)
})
context("test-plot_regional_similarity")



# Load local_cossim object
regional_sims <- readRDS(system.file("states/regional_sims.rds",
  package = "MutationalPatterns"
))

# Plot the regional similarity
output = plot_regional_similarity(regional_sims)

# Plot outlier samples with a different color.
output_outlier = plot_regional_similarity(regional_sims, max_cossim = 0.5)

# Plot samples per chromosome
output_l = plot_regional_similarity(regional_sims, per_chrom = TRUE)

# Run tests
test_that("Output has correct class", {
    expect_true(inherits(output, c("gg")))
    expect_true(inherits(output_outlier, c("gg")))
    expect_true(inherits(output_l, c("list")))
    expect_true(inherits(output_l[[1]], c("gg")))
})

test_that("Output per chromosome has correct length", {
    expect_equal(length(output_l), 6)
})
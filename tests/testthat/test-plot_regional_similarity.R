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

# Plot samples with a rug
output_rug = plot_regional_similarity(regional_sims, plot_rug = TRUE)

# Use custom x-axis breaks
output_xbreaks = plot_regional_similarity(regional_sims, x_axis_breaks = c(30, 66, 300))

# Run tests
test_that("Output has correct class", {
    expect_true(inherits(output, c("gg")))
    expect_true(inherits(output_outlier, c("gg")))
    expect_true(inherits(output_l, c("list")))
    expect_true(inherits(output_l[[1]], c("gg")))
    expect_true(inherits(output_rug, c("gg")))
    expect_true(inherits(output_xbreaks, c("gg")))
})

test_that("Output per chromosome has correct length", {
    expect_equal(length(output_l), 3)
})
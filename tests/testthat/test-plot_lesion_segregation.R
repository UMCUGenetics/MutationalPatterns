context("test-plot_lesion_segregation")


# Load GRangesList
grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

# Only use two samples to reduce runtime
grl <- grl[1:2]

# Select sample
gr <- grl[[1]]

# Perform function
output <- plot_lesion_segregation(grl)
output_singlesample <- plot_lesion_segregation(gr, sample_name = "Colon1")
output_noname <- plot_lesion_segregation(gr)
output_l <- plot_lesion_segregation(gr, per_chrom = TRUE, sample_name = "Colon1")
output_chr_filter = plot_lesion_segregation(grl, chromosomes = c("chr2", "chr3"))
output_chr_filter2 = plot_lesion_segregation(grl, chromosomes = c("2", "3"))
output_subsample <- plot_lesion_segregation(grl, subsample = 0.1)

test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_singlesample, c("gg")))
  expect_true(inherits(output_noname, c("gg")))
  expect_true(inherits(output_l, c("list")))
  expect_true(inherits(output_l[[1]], c("gg")))
  expect_true(inherits(output_chr_filter, c("gg")))
  expect_true(inherits(output_chr_filter2, c("gg")))
  expect_true(inherits(output_subsample, c("gg")))
})

test_that("Output per chromosome has correct length", {
  expect_equal(length(output_l), 22)
})

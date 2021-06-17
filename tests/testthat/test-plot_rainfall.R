context("test-plot_rainfall")

# Read data
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

grl_indel_context <- readRDS(system.file("states/blood_grl_indel_context.rds",
  package = "MutationalPatterns"
))

grl_dbs_context <- readRDS(system.file("states/blood_grl_dbs_context.rds",
  package = "MutationalPatterns"
))

mbs_grl <- readRDS(system.file("states/blood_grl_mbs.rds",
  package = "MutationalPatterns"
))


# Specify chromosomes of interest.
chromosomes <- paste0("chr", c(1:22))

# Do a rainfall plot for all chromosomes:
output <- plot_rainfall(vcfs[[1]], title = names(vcfs[1]), chromosomes = chromosomes)

# Plot a single chromosome (chromosome 1):
output_singlechrom <- plot_rainfall(vcfs[[1]], title = names(vcfs[1]), chromosomes = chromosomes[1])

# Plot a subset of the variants
output_subset <- plot_rainfall(vcfs[[1]][1:10], title = names(vcfs[1]), chromosomes = chromosomes)

# plot indels
output_indel <- plot_rainfall(grl_indel_context[[1]], title = names(grl_indel_context[1]), chromosomes = chromosomes, type = "indel")

# plot dbs
output_dbs <- plot_rainfall(grl_dbs_context[[1]], title = names(grl_dbs_context[1]), chromosomes = chromosomes, type = "dbs")

# plot mbs
output_mbs <- plot_rainfall(mbs_grl[[1]], title = names(mbs_grl[1]), chromosomes = chromosomes, type = "mbs")

# Plot an empty gr
output_empty <- plot_rainfall(vcfs[[1]][0], title = names(vcfs[1]), chromosomes = chromosomes)


test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_singlechrom, c("gg")))
  expect_true(inherits(output_subset, c("gg")))
  expect_true(inherits(output_indel, c("gg")))
  expect_true(inherits(output_dbs, c("gg")))
  expect_true(inherits(output_mbs, c("gg")))
  expect_true(inherits(output_empty, c("gg")))
})

test_that("Subsetted output contains the correct subset of colours", {
  colours_used <- unique(ggplot_build(output_subset)$data[[1]][["colour"]])
  expect_equal(colours_used, c("#ADCC54", "#DE1C14", "#2EBAED", "#D4D2D2"))
})

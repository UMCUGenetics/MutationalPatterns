context("test-determine_regional_similarity")

# See the 'read_vcfs_as_granges()' example for how we obtained the
# following data:
grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

# We pool all the variants together, because the function doesn't work well with a limited number of mutations.
# Still, in practice we recommend to use more mutations that in this example.
gr = unlist(grl)

# Specifiy the chromosomes of interest.
chromosomes <- names(genome(gr)[1:6])

# Load the corresponding reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Determine the regional similarities. Here we use a small window size to make the function work.
# In practice, we recommend a larger window size.
output = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 40, stepsize = 10, max_window_size_gen = 40000000)
output_tri = determine_regional_similarity(gr, ref_genome, chromosomes[6], window_size = 40, stepsize = 40, tri_correction = TRUE, max_window_size_gen = 40000000)
output_notexcl = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 40, stepsize = 10, exclude_self_mut_mat = FALSE, max_window_size_gen = 40000000)

# Load expected
expected <- readRDS(system.file("states/regional_sims.rds",
  package = "MutationalPatterns"
))

# Run tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("region_cossim")))
  expect_true(inherits(output_tri, c("region_cossim")))
  expect_true(inherits(output_notexcl, c("region_cossim")))
  
})

test_that("Output has correct dimensions", {
  expect_equal(dim(output@sim_tb), c(152, 8))
  expect_equal(dim(output@pos_tb), c(1789, 3))
  expect_equal(dim(output_tri@sim_tb), c(6, 10))
  expect_equal(dim(output_tri@pos_tb), c(269, 3))
  expect_equal(dim(output_notexcl@sim_tb), c(152, 8))
  expect_equal(dim(output_notexcl@pos_tb), c(1789, 3))
})

test_that("transforms correctly", {
  expect_equal(output, expected)
})

test_that("exclude_self_mut_mat reduces cosine similarity", {
  expect_equal(sum(output@sim_tb$cossim > output_notexcl@sim_tb$cossim), 0)
})

test_that("oligonucleotide frequency correction increases cosine similarity", {
  expect_equal(sum(output_tri@sim_tb$corrected_cossim < output_tri@sim_tb$cossim), 0)
})

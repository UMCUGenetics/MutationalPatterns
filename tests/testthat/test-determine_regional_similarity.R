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
chromosomes <- names(genome(gr)[1:3])

# Load the corresponding reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Determine the regional similarities. Here we use a small window size to make the function work.
# In practice, we recommend a larger window size.
output = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 40, stepsize = 10, max_window_size_gen = 40000000)
output_oligo = determine_regional_similarity(gr, ref_genome, chromosomes[3], window_size = 60, stepsize = 60, oligo_correction = TRUE, max_window_size_gen = 40000000)
output_notexcl = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 40, stepsize = 10, exclude_self_mut_mat = FALSE, max_window_size_gen = 40000000)
output_smallext = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 40, stepsize = 10, extension = 0, max_window_size_gen = 40000000)
output_verbose = determine_regional_similarity(gr, ref_genome, chromosomes, window_size = 40, stepsize = 40, max_window_size_gen = 40000000, verbose = TRUE)

# Load expected
expected <- readRDS(system.file("states/regional_sims.rds",
                                package = "MutationalPatterns"
))

# Run tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("region_cossim")))
  expect_true(inherits(output_oligo, c("region_cossim")))
  expect_true(inherits(output_notexcl, c("region_cossim")))
  expect_true(inherits(output_smallext, c("region_cossim")))
  expect_true(inherits(output_verbose, c("region_cossim")))
})

test_that("Output has correct dimensions", {
  expect_equal(dim(output@sim_tb), c(76, 8))
  expect_equal(dim(output@pos_tb), c(915, 3))
  expect_equal(dim(output_oligo@sim_tb), c(2, 10))
  expect_equal(dim(output_oligo@pos_tb), c(294, 3))
  expect_equal(dim(output_notexcl@sim_tb), c(76, 8))
  expect_equal(dim(output_notexcl@pos_tb), c(915, 3))
  expect_equal(dim(output_smallext@sim_tb), c(76, 8))
  expect_equal(dim(output_smallext@pos_tb), c(915, 3))
  expect_equal(dim(output_verbose@sim_tb), c(21, 8))
  expect_equal(dim(output_verbose@pos_tb), c(915, 3))
})

test_that("transforms correctly", {
  expect_equal(output@sim_tb, expected@sim_tb, check.attributes = FALSE)
  expect_equal(output@pos_tb, expected@pos_tb, check.attributes = FALSE)
})

test_that("exclude_self_mut_mat reduces cosine similarity", {
  expect_equal(sum(output@sim_tb$cossim > output_notexcl@sim_tb$cossim), 0)
})

test_that("oligonucleotide frequency correction increases cosine similarity", {
  expect_equal(sum(output_oligo@sim_tb$corrected_cossim < output_oligo@sim_tb$cossim), 0)
})

test_that("A smaller extension increases cosine similarity", {
  expect_equal(sum(output@sim_tb$cossim > output_smallext@sim_tb$cossim), 0)
})

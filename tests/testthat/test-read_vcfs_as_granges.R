context("test-read_vcfs_as_granges")

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome)
library(ref_genome, character.only = TRUE)

sample_names <- c(
  "colon1", "colon2", "colon3",
  "intestine1", "intestine2", "intestine3",
  "liver1", "liver2", "liver3"
)

vcfs <- list.files(system.file("extdata", package = "MutationalPatterns"),
  pattern = "sample.vcf", full.names = TRUE
)

# Test default
test_that("loads multiple samples", {
  output <- read_vcfs_as_granges(vcfs, sample_names, ref_genome)
  expect_that(length(output), equals(9))
  expect_true(inherits(output, "CompressedGRangesList"))
})

# Test for seqlevel filters
test_that("nuclear filter works", {
  output <- read_vcfs_as_granges(vcfs, sample_names, ref_genome)
  expected <- c(
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
    "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY"
  )

  expect_that(seqlevels(output), equals(expected))
})

test_that("autosomal filter works", {
  output <- read_vcfs_as_granges(vcfs, sample_names, ref_genome, "auto")
  expected <- c(
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
    "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22"
  )

  expect_that(seqlevels(output), equals(expected))
})

test_that("unfiltered works", {

  ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  library(ref_genome, character.only = TRUE)

  output <- read_vcfs_as_granges(vcfs, sample_names, ref_genome, "none")
  expected <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
                "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
                "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
                "chrX", "chrY", "chrM", "GL000207.1", "GL000226.1", "GL000229.1", 
                "GL000231.1", "GL000210.1", "GL000239.1", "GL000235.1", "GL000201.1", 
                "GL000247.1", "GL000245.1", "GL000197.1", "GL000203.1", "GL000246.1", 
                "GL000249.1", "GL000196.1", "GL000248.1", "GL000244.1", "GL000238.1", 
                "GL000202.1", "GL000234.1", "GL000232.1", "GL000206.1", "GL000240.1", 
                "GL000236.1", "GL000241.1", "GL000243.1", "GL000242.1", "GL000230.1", 
                "GL000237.1", "GL000233.1", "GL000204.1", "GL000198.1", "GL000208.1", 
                "GL000191.1", "GL000227.1", "GL000228.1", "GL000214.1", "GL000221.1", 
                "GL000209.1", "GL000218.1", "GL000220.1", "GL000213.1", "GL000211.1", 
                "GL000199.1", "GL000217.1", "GL000216.1", "GL000215.1", "GL000205.1", 
                "GL000219.1", "GL000224.1", "GL000223.1", "GL000195.1", "GL000212.1", 
                "GL000222.1", "GL000200.1", "GL000193.1", "GL000194.1", "GL000225.1", 
                "GL000192.1")
  
  expect_equal(seqlevels(output), expected)
})


# Test that a warning is given when vcf and names lengths don't match
test_that("An error is given when vcf and names lengths don't match", {
  expect_error(
    {
      read_vcfs_as_granges(vcfs, sample_names[1:8], ref_genome)
    },
    "Please provide the same number of sample names as VCF files"
  )
})

# Test that a warning is given when the supplied reference is not a BSgenome object
test_that("An error is given when the supplied ref is not a BSgenome object", {
  expect_error(
    {
      read_vcfs_as_granges(vcfs, sample_names, "a")
    },
    "Please provide the name of a BSgenome object."
  )
})

# Test that you can read in specific mutation types
vcf_fnames <- list.files(system.file("extdata", package = "MutationalPatterns"),
  pattern = "blood.*vcf", full.names = TRUE
)
sample_names <- c("AC", "ACC55", "BCH")
test_that("indels work", {
  output <- read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "indel")
  expect_that(length(output), equals(3))
  expect_true(inherits(output, "CompressedGRangesList"))
})

test_that("dbs work", {
  output <- read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "dbs")
  expect_that(length(output), equals(3))
  expect_true(inherits(output, "CompressedGRangesList"))
})

output <- read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "mbs")
test_that("mbs work", {
  expect_that(length(output), equals(3))
  expect_true(inherits(output, "CompressedGRangesList"))
})

test_that("all mutation types work", {
  output <- read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "all")
  expect_that(length(output), equals(3))
  expect_true(inherits(output, "CompressedGRangesList"))
})

test_that("predefined_dbs_mbs argument works", {
  output <- read_vcfs_as_granges(vcf_fnames, 
                                 sample_names, 
                                 ref_genome, 
                                 type = "dbs",
                                 predefined_dbs_mbs = TRUE)
  expect_that(length(output), equals(3))
  expect_that(as.double(S4Vectors::elementNROWS(output)), equals(c(0, 0, 0)))
  expect_true(inherits(output, "CompressedGRangesList"))
})

# Test function works on an empty vcf
empty_vcf <- list.files(system.file("extdata", package = "MutationalPatterns"),
  pattern = "empty.vcf", full.names = TRUE
)

test_that("Empty vcf works", {
  expect_warning(
    {
      output <- read_vcfs_as_granges(empty_vcf, "empty", ref_genome)
    },
    "There were 0 variants \\(before filtering\\) found in the vcf file"
  )
  expect_true(inherits(output, "CompressedGRangesList"))
  expect_equal(length(output[[1]]), 0)
})

# Test that it is possible to keep duplicate variants
duplicate_vcf <- list.files(system.file("extdata", package = "MutationalPatterns"),
                        pattern = "duplicate.vcf", full.names = TRUE
)

test_that("Empty vcf works", {
  output <- read_vcfs_as_granges(duplicate_vcf, "duplicate", ref_genome, remove_duplicate_variants = F)
  expect_true(inherits(output, "CompressedGRangesList"))
  expect_equal(length(output[[1]]), 2)
})

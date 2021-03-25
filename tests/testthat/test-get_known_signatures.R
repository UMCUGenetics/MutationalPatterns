context("test-get_known_signatures")

# Get reference snv signature from COSMIC
output <- get_known_signatures()

# Get reference snv signature from COSMIC,
# including potential artifacts.
output_artifact <- get_known_signatures(incl_poss_artifacts = TRUE)

# Get a GRCh38 version of the signatures
output_grch38 = get_known_signatures(genome = "GRCh38")

# Get dbs signatures
output_dbs <- get_known_signatures("dbs")

# Get indel signatures
output_indel <- get_known_signatures("indel")

# Get transcription strand bias snv signatures
output_tsb <- get_known_signatures("tsb_snv")

# Get a GRCh38 version of the signatures
output_v3_1 = get_known_signatures(source = "COSMIC_v3.1")

# Get reference signatures from SIGNAL
output_signal <- get_known_signatures(source = "SIGNAL")

# Get reference signatures from SIGNAL,
# including potential artifacts
output_signal_arti <- get_known_signatures(source = "SIGNAL", incl_poss_artifacts = TRUE)

## Get exposure signatures from SIGNAL
output_exposure <- get_known_signatures(source = "SIGNAL", sig_type = "exposure")

# Get DBS exposure signatures from SIGNAL
output_exposure_dbs <- get_known_signatures("dbs", source = "SIGNAL", sig_type = "exposure")

# Get all tissue specific signatures from SIGNAL
output_tissue <- get_known_signatures(source = "SIGNAL", sig_type = "tissue")


# Get Bladder specific signatures from SIGNAL
output_bladder <- get_known_signatures(
  source = "SIGNAL",
  sig_type = "tissue",
  tissue_type = "Bladder"
)

# Get sparse signatures
output_sparse <- get_known_signatures(source = "SPARSE")


# Perform tests
test_that("Output has correct class", {
  expect_true(inherits(output, "matrix"))
  expect_true(inherits(output_artifact, "matrix"))
  expect_true(inherits(output_grch38, "matrix"))
  expect_true(inherits(output_dbs, "matrix"))
  expect_true(inherits(output_indel, "matrix"))
  expect_true(inherits(output_tsb, "matrix"))
  expect_true(inherits(output_v3_1, "matrix"))
  expect_true(inherits(output_signal, "matrix"))
  expect_true(inherits(output_signal_arti, "matrix"))
  expect_true(inherits(output_exposure, "matrix"))
  expect_true(inherits(output_exposure_dbs, "matrix"))
  expect_true(inherits(output_tissue, "matrix"))
  expect_true(inherits(output_bladder, "matrix"))
  expect_true(inherits(output_sparse, "matrix"))
})

test_that("Output has correct dimensions", {
  expect_equal(dim(output), c(96, 60))
  expect_equal(dim(output_artifact), c(96, 78))
  expect_equal(dim(output_grch38), c(96, 60))
  expect_equal(dim(output_dbs), c(78, 11))
  expect_equal(dim(output_indel), c(83, 18))
  expect_equal(dim(output_tsb), c(192, 65))
  expect_equal(dim(output_v3_1), c(96, 54))
  expect_equal(dim(output_signal), c(96, 29))
  expect_equal(dim(output_signal_arti), c(96, 41))
  expect_equal(dim(output_exposure), c(96, 54))
  expect_equal(dim(output_exposure_dbs), c(78, 8))
  expect_equal(dim(output_tissue), c(96, 192))
  expect_equal(dim(output_bladder), c(96, 7))
  expect_equal(dim(output_sparse), c(96, 10))
})


test_that("An error is thrown when tissue_type is used outside tissue mode.", {
  expect_error(
    {
      get_known_signatures(tissue_type = "Bladder")
    },
    "tissue_type can only be used with `sig_type = 'tissue'`"
  )
})

test_that("An error is thrown when a file doesn't exist.", {
  expect_error(
    {
      get_known_signatures("dbs", source = "SIGNAL")
    },
    paste0(
      "Look at the documentation of \\'get_known_signatures\\(\\)\\' for",
      " all the possible combinations of arguments."
    )
  )
})

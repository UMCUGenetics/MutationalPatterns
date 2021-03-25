context("test-merge_signatures")


# Get signatures
signatures <- get_known_signatures()

## Merge signatures
output <- merge_signatures(signatures)


## Merge signatures using a stricter cutoff
output_strict <- merge_signatures(signatures, cos_sim_cutoff = 0.9)

## Merge signatures using a different merging character
output_char <- merge_signatures(signatures, merge_char = "_")

# Tests
test_that("Output has correct class", {
  expect_true(inherits(output, "matrix"))
  expect_true(inherits(output_strict, "matrix"))
  expect_true(inherits(output_char, "matrix"))
})


test_that("Output sizes are correct", {
  expect_equal(dim(output), c(96, 45))
  expect_equal(dim(output_strict), c(96, 58))
  expect_equal(dim(output_char), c(96, 45))
})

test_that("Signature names are correct", {
  expect_equal(colnames(output), c("SBS1", "SBS2", "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS9", 
                                   "SBS10b", "SBS11", "SBS13", "SBS14", "SBS16", "SBS17a", "SBS17b", 
                                   "SBS20", "SBS21", "SBS22", "SBS25", "SBS28", "SBS30", "SBS32", 
                                   "SBS33", "SBS34", "SBS35", "SBS38", "SBS41", "SBS42", "SBS44", 
                                   "SBS84", "SBS85", "SBS87", "SBS88", "SBS89", "SBS90", "SBS91", 
                                   "SBS93", "SBS36;SBS18", "SBS10d;SBS10a;SBS10c", "SBS15;SBS6", 
                                   "SBS29;SBS24", "SBS26;SBS12;SBS37", "SBS86;SBS39", "SBS40;SBS3;SBS92;SBS5", 
                                   "SBS94;SBS4;SBS8", "SBS23;SBS19;SBS31"))
  expect_equal(colnames(output_strict), c("SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS6", "SBS7a", "SBS7b", 
                                          "SBS7c", "SBS7d", "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
                                          "SBS10d", "SBS11", "SBS13", "SBS14", "SBS15", "SBS16", "SBS17a", 
                                          "SBS17b", "SBS19", "SBS20", "SBS21", "SBS22", "SBS23", "SBS24", 
                                          "SBS25", "SBS28", "SBS29", "SBS30", "SBS31", "SBS32", "SBS33", 
                                          "SBS34", "SBS35", "SBS37", "SBS38", "SBS39", "SBS40", "SBS41", 
                                          "SBS42", "SBS44", "SBS84", "SBS85", "SBS86", "SBS87", "SBS88", 
                                          "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", "SBS94", "SBS26;SBS12", 
                                          "SBS36;SBS18"))
  expect_equal(colnames(output_char), c("SBS1", "SBS2", "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS9", 
                                        "SBS10b", "SBS11", "SBS13", "SBS14", "SBS16", "SBS17a", "SBS17b", 
                                        "SBS20", "SBS21", "SBS22", "SBS25", "SBS28", "SBS30", "SBS32", 
                                        "SBS33", "SBS34", "SBS35", "SBS38", "SBS41", "SBS42", "SBS44", 
                                        "SBS84", "SBS85", "SBS87", "SBS88", "SBS89", "SBS90", "SBS91", 
                                        "SBS93", "SBS36_SBS18", "SBS10d_SBS10a_SBS10c", "SBS15_SBS6", 
                                        "SBS29_SBS24", "SBS26_SBS12_SBS37", "SBS86_SBS39", "SBS40_SBS3_SBS92_SBS5", 
                                        "SBS94_SBS4_SBS8", "SBS23_SBS19_SBS31"))
})

# Test that an error is given when the signature names already contain the merge_char
signatures_colnames <- signatures
colnames(signatures_colnames)[1] <- "SBS;1"
test_that("An error is given when signature names contain merge_char", {
  expect_error(
    {
      merge_signatures(signatures_colnames)
    },
    "Please remove all ; characters from your signature names."
  )
})

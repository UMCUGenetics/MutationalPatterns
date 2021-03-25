library(tidyverse)
library(VariantAnnotation)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Function to add fake MNVs to example data. This is usefull for the unit tests, examples and vignette.
add_fake_mnvs <- function(vcf, ref_genome, nr_mnvs = NA) {
  if (is.na(nr_mnvs)) {
    nr_mnvs <- sample(seq(5, 10), 1)
  }

  mnv_lengths <- sample(seq(3, 10),
    nr_mnvs,
    replace = TRUE,
    prob = c(0.35, 0.25, 0.15, 0.1, 0.05, 0.05, 0.025, 0.025)
  )
  mnvs <- purrr::map(mnv_lengths, create_1_mnv, vcf) %>%
    do.call(rbind, .) %>%
    sort()
  vcf <- rbind(vcf, mnvs)
  return(vcf)
}

# helper function of add_fake_mnvs
create_1_mnv <- function(mnv_length, vcf) {

  # Select random variant
  vcf <- vcf[isSNV(vcf)]
  variant <- vcf[sample.int(length(vcf), 1)]

  # Combine variant `mnv_length` times
  mnv_variant <- purrr::map(seq(1, mnv_length), function(x) {
      return(variant)
    }) %>%
    do.call(rbind, .)

  # Replace the positions in the vcf
  first_start <- start(mnv_variant)[1] + 100
  starts <- seq(first_start, first_start + (mnv_length - 1))
  end(mnv_variant) <- starts
  start(mnv_variant) <- starts

  # Update the ref bases
  gr <- granges(mnv_variant)
  old_seqnames <- as.vector(seqnames(gr))
  seqlevelsStyle(gr) <- "UCSC"
  seqlevels(gr) <- str_c("chr", c(1:22, "X", "Y"))
  refs <- Biostrings::getSeq(BSgenome::getBSgenome(ref_genome), gr)
  ref(mnv_variant) <- refs

  possible_alts <- purrr::map(as.vector(refs), ~ setdiff(c("A", "C", "T", "G"), .x))
  alts <- purrr::map(possible_alts, ~ sample(.x, 1)) %>%
    DNAStringSetList(use.names = FALSE)
  alt(mnv_variant) <- alts

  # Update the names of the new variants
  names(mnv_variant) <- str_c(
    old_seqnames,
    ":",
    starts,
    "_",
    as.vector(ref(mnv_variant)),
    "/",
    as.vector(unlist(alt(mnv_variant)))
  )
  return(mnv_variant)
}


# Function to remove data from a vcf, that is not needed for the package.
decrease_vcf_size <- function(vcf_fname) {
  # Read vcf
  vcf <- readVcf(vcf_fname)

  # Remove unnecessary info
  info(vcf) <- info(vcf)[, 0, drop = FALSE]
  info(header(vcf)) <- info(header(vcf))[0, , drop = FALSE]

  # Remove unnecessary genotype data
  gt_to_keep <- c("GT", "AD", "DP", "GQ", "PL")
  geno(vcf) <- geno(vcf)[gt_to_keep]
  geno(header(vcf)) <- geno(header(vcf))[gt_to_keep, , drop = FALSE]

  # Remove filter statements from header
  fixed(header(vcf))[["FILTER"]] <- fixed(header(vcf))[["FILTER"]][0, , drop = FALSE]

  # Remove all but one sample
  tmp_name <- "tmp.vcf"
  tmp_name2 <- "tmp2.vcf"
  writeVcf(vcf, tmp_name)
  system(str_c("grep -v '##ALT' ", tmp_name, " > ", tmp_name2))
  system(str_c("cut -f1-10 ", tmp_name2, " > ", vcf_fname))
  file.remove(tmp_name)
  file.remove(tmp_name2)
  invisible(0)
}


# Determine vcf and donor names
snv_vcf_fnames <- list.files("~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/SNVs/hg19/Healthy_bone_marrow/",
  pattern = "MQ60.vcf",
  full.names = TRUE
)
indel_vcf_fnames <- list.files("~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/INDELs/hg19/Healthy_bone_marrow/",
  pattern = "CALLABLE.vcf",
  full.names = TRUE
)
vcf_fnames <- c(snv_vcf_fnames, indel_vcf_fnames)
sample_names <- vcf_fnames %>%
  basename() %>%
  str_remove("_.*")
donor_names <- sample_names %>%
  str_remove("HSC.*") %>%
  str_remove("MPP.*")

used_donors <- c("AC", "ACC55", "BCH")
vcf_f <- donor_names %in% used_donors
vcf_fnames <- vcf_fnames[vcf_f]
donor_names <- donor_names[vcf_f]

# Read vcfs and combine vcfs from the same donor
vcf_l <- purrr::map(vcf_fnames, readVcf, genome = "hg19")
vcf_l <- vcf_l %>%
  split(donor_names) %>% # Split per donor into list of lists
  purrr::map(., function(x) do.call(rbind, x)) %>% # Combine vcfs from a single donor
  purrr::map(sort) %>%
  purrr::map(unique)

# Add fake mnvs.
set.seed(42)
vcf_l <- purrr::map(vcf_l, add_fake_mnvs, ref_genome)

# Write vcfs
out_vcf_fnames <- str_c("inst/extdata/blood-", names(vcf_l), ".vcf")
purrr::map2(vcf_l, out_vcf_fnames, writeVcf)

# Decrease vcf size
purrr::map(out_vcf_fnames, decrease_vcf_size)

# Create granges object from the vcfs
vcf_fnames <- list.files("inst/extdata/", pattern = "blood.*.vcf", full.names = TRUE)
vcf_l <- purrr::map(vcf_fnames, readVcf, "hg19")
sample_names <- basename(vcf_fnames) %>%
  str_remove("blood-") %>%
  str_remove(".vcf")
names(vcf_l) <- sample_names
grl <- purrr::map(vcf_l, granges) %>%
  GRangesList()
seqlevelsStyle(grl) <- "UCSC"
seqlevels(grl, pruning.mode = "fine") <- str_c("chr", c(1:22, "X"))
saveRDS(grl, "inst/states/blood_grl.rds")


# Create empty vcf
vcf <- vcf_l[[1]]
vcf <- vcf[0]
writeVcf(vcf, "inst/extdata/empty.vcf")


# Reduce size of old example vcfs
samples = c("colon1", "colon2", "colon3", "intestine1", "intestine2", "intestine3", "liver1", "liver2", "liver3")
vcf_fnames = paste0("inst/extdata/", samples, "-sample.vcf")
purrr::map(vcf_fnames, decrease_vcf_size)

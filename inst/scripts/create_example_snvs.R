library(tidyverse)
library(VariantAnnotation)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

vcf_files <- list.files(system.file("extdata", package = "MutationalPatterns"),
  pattern = "sample.vcf", full.names = TRUE
)

sample_names <- c(
  "colon1", "colon2", "colon3",
  "intestine1", "intestine2", "intestine3",
  "liver1", "liver2", "liver3"
)

# Create grl
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
saveRDS(grl, "inst/states/read_vcfs_as_granges_output.rds")

# Create snvs
mut_mat <- mut_matrix(grl, ref_genome)
saveRDS(mut_mat, "inst/states/mut_mat_data.rds")

mut_mat_extended <- mut_matrix(grl, ref_genome, extension = 2)
saveRDS(mut_mat_extended, "inst/states/mut_mat_data_extended.rds")


# Create  transcription strand matrix
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg19)
saveRDS(mut_mat_s, "inst/states/mut_mat_s_data.rds")

# Create replication direction
repli_file <- system.file("extdata/ReplicationDirectionRegions.bed",
  package = "MutationalPatterns"
)
repli_strand <- read.table(repli_file, header = TRUE)
repli_strand_granges <- GRanges(
  seqnames = repli_strand$Chr,
  ranges = IRanges(
    start = repli_strand$Start + 1,
    end = repli_strand$Stop
  ),
  strand_info = repli_strand$Class
)
seqlevelsStyle(repli_strand_granges) <- "UCSC"
saveRDS(repli_strand_granges, "inst/states/repli_strand.rds")

# Create replication strand matrix
mut_mat_repli <- mut_matrix_stranded(grl, ref_genome, repli_strand_granges, mode = "replication")
saveRDS(mut_mat_repli, "inst/states/mut_mat_repli.rds")

# Extract signatures
nmf_res <- extract_signatures(mut_mat, rank = 2)
saveRDS(nmf_res, "inst/states/nmf_res_data.rds")
nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2)
saveRDS(nmf_res_strand, "inst/states/nmf_res_strand_data.rds")

# Get signatures
signatures <- get_known_signatures()

# Normal refit
fit_res <- fit_to_signatures(mut_mat, signatures)
saveRDS(fit_res, "inst/states/snv_refit.rds")

# Strict refit
strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.05)
saveRDS(strict_refit$fit_res, "inst/states/strict_snv_refit.rds")

strict_refit_best <- fit_to_signatures_strict(mut_mat, signatures[,1:5], max_delta = 0.004, method = "best_subset")
saveRDS(strict_refit_best$fit_res, "inst/states/strict_best_snv_refit.rds")


# bootstrapped refit
set.seed(42)
contri_boots <- fit_to_signatures_bootstrapped(mut_mat, signatures, n_boots = 2, max_delta = 0.05)
saveRDS(contri_boots, "inst/states/bootstrapped_snv_refit.rds")

# Calculate lesion segregation
lesion_segretation <- calculate_lesion_segregation(grl[1:2], sample_names[1:2])
saveRDS(lesion_segretation, "inst/states/lesion_segregation.rds")


# Split mutation types
# Read in genomic regions
CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
  package = "MutationalPatterns"
))
promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
  package = "MutationalPatterns"
))
flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
  package = "MutationalPatterns"
))

# Combine the regions into a single GRangesList
regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
names(regions) <- c("Promoter", "Promoter flanking", "CTCF")

seqlevelsStyle(regions) <- "UCSC"

# Read in some variants.
grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

grl_split <- split_muts_region(grl, regions)
saveRDS(grl_split, "inst/states/grl_split_region.rds")

mut_mat_split_region <- mut_matrix(grl_split, ref_genome)
saveRDS(mut_mat_split_region, "inst/states/mut_mat_splitregions.rds")

mut_mat_longregion <- lengthen_mut_matrix(mut_mat_split_region)
saveRDS(mut_mat_longregion, "inst/states/mut_mat_longregions.rds")



# Create context potential damage tibble
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

contexts = rownames(mut_mat)[1:6]
gene_ids = c(7157)
context_mismatches = context_potential_damage_analysis(contexts, txdb, ref_genome, gene_ids)
saveRDS(context_mismatches, "inst/states/context_mismatches.rds")


# Specifiy the chromosomes of interest.
chromosomes <- names(genome(grl)[1:3])
regional_sims = determine_regional_similarity(unlist(grl), ref_genome, chromosomes, window_size = 40, stepsize = 10, max_window_size_gen = 40000000)
saveRDS(regional_sims, "inst/states/regional_sims.rds")


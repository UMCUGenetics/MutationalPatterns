#' Checks that there are no variants with multiple alt alleles.
#'
#' @details
#' This function checks for a GRanges/GRangesList object, that there are no variants with multiple alternative alleles.
#' It works by checking that the number of variants is equal to the number of alternative alleles.
#' It trows an error when this is not the case.
#'
#' @param grl GRanges/GRangesList object
#'
#' @return Invisibly returns its input when it succeeds.
#' @noRd
#' @importFrom magrittr %>%

.check_no_multi_alts <- function(grl) {
  if (inherits(grl, "CompressedGRangesList")) {
    gr <- unlist(grl)
  } else if (inherits(grl, "GRanges")) {
    gr <- grl
  } else {
    .not_gr_or_grl(grl)
  }
  .check_no_multi_alts_gr(gr)
  invisible(grl)
}

#' Checks that there are no variants with multiple alt alleles.
#'
#' @details
#' This function checks for a single GRanges object, that there are no variants with multiple alternative alleles.
#' It works by checking that the number of variants is equal to the number of alternative alleles.
#' It trows an error when this is not the case.
#'
#' @param gr GRanges object
#'
#' @return Invisibly returns its input when it succeeds.
#' @noRd
#' @importFrom magrittr %>%

.check_no_multi_alts_gr <- function(gr) {
  alt <- .get_alt(gr)
  nr_alts <- alt %>%
    unlist() %>%
    length()
  if (length(gr) != nr_alts) {
    stop(paste0(
      "There should not be any variants with multiple alternative alleles.\n",
      "You can remove these by using the `get_mut_type` function."
    ), call. = FALSE)
  }
  invisible(gr)
}

#' Checks that there are no SNVs/MNVs.
#'
#' @details
#' This function checks for a GRanges/GRangesList object, that there are no SNVs/MNVs.
#' It trows an error when this is not the case.
#'
#' @param grl GRanges/GrangesList object
#' @noRd
#' @return Invisibly returns its input when it succeeds.
#'

.check_no_substitutions <- function(grl) {
  if (inherits(grl, "CompressedGRangesList")) {
    gr <- unlist(grl)
  } else if (inherits(grl, "GRanges")) {
    gr <- grl
  } else {
    .not_gr_or_grl(grl)
  }
  .check_no_substitutions_gr(gr)
  invisible(grl)
}

#' Checks that there are no SNVs/MNVs.
#'
#' @details
#' This function checks for a single GRanges object, that there are no SNVs/MNVs.
#' It trows an error when this is not the case.
#'
#' @param gr GRanges object
#' @noRd
#' @return Invisibly returns its input when it succeeds.
#'

.check_no_substitutions_gr <- function(gr) {
  snv_f <- .find_substitution(gr)
  nr_snv <- sum(snv_f)
  snv_present <- nr_snv >= 1
  if (snv_present) {
    stop(stringr::str_c("There seem to be ", 
                        nr_snv, 
                        " substitutions present in your data. ",
                        "Make sure to remove all substitutions with the `get_mut_type` function."), 
         call. = FALSE)
  }
  invisible(gr)
}


#' Checks that there are no SNVs/MNVs.
#'
#' @details
#' This function checks for a GRanges/GRangesList object, that there are no Indels.
#' It trows an error when this is not the case.
#'
#' @param grl GRanges/GrangesList object
#' @noRd
#' @return Invisibly returns its input when it succeeds.
#'
.check_no_indels <- function(grl) {
  if (inherits(grl, "CompressedGRangesList")) {
    gr <- unlist(grl)
  } else if (inherits(grl, "GRanges")) {
    gr <- grl
  } else {
    .not_gr_or_grl(grl)
  }
  .check_no_indels_gr(gr)
  invisible(grl)
}

#' Checks that there are no Indels.
#'
#' @details
#' This function checks for a single GRanges object, that there are no Indels.
#' It trows an error when this is not the case.
#'
#' @param gr GRanges object
#' @noRd
#' @return Invisibly returns its input when it succeeds.
#'
.check_no_indels_gr <- function(gr) {
  snv_f <- .find_substitution(gr)
  nr_indel <- sum(!snv_f)
  snv_present <- nr_indel >= 1
  if (snv_present) {
    stop(paste0(
      "There seem to be at least ", nr_indel,
      " Indels present in your data. ",
      "Make sure to remove all Indels with the `get_mut_type` function."
    ), call. = FALSE)
  }
  invisible(gr)
}


#' Removes variants with multiple alt alleles.
#'
#' @details
#' This function removes variants with multiple alternative alleles for a GRanges/GRangesList object.
#'
#' @param grl GRanges/GRangesList object
#' @noRd
#' @return A filtered version of the input GRanges/GRangesList object.
.remove_multi_alts_variants <- function(grl) {
  if (inherits(grl, "CompressedGRangesList")) {
    gr_l <- as.list(grl)
    grl <- purrr::map(gr_l, .remove_multi_alts_variants_gr) %>%
      GRangesList()
    return(grl)
  } else if (inherits(grl, "GRanges")) {
    gr <- .remove_multi_alts_variants_gr(grl)
    return(gr)
  } else {
    .not_gr_or_grl(grl)
  }
}

#' Removes variants with multiple alt alleles.
#'
#' @details
#' This function removes variants with multiple alternative alleles for a single GRanges object.
#'
#' @param gr GRanges object
#' @noRd
#' @return A filtered version of the input GRanges object.
#'
.remove_multi_alts_variants_gr <- function(gr) {
  alt <- .get_alt(gr)
  gr <- gr[IRanges::elementNROWS(alt) == 1]
  return(gr)
}

#' Removes SNV/MNV variants.
#'
#' @details
#' This function removes SNV/MNV variants for a GRanges/GrangesList object.
#'
#' @param grl GRanges/GRangesList object
#' @noRd
#' @return A filtered version of the input GRanges/GrangesList object.
.remove_substitutions <- function(grl) {
  if (inherits(grl, "CompressedGRangesList")) {
    gr_l <- as.list(grl)
    grl <- purrr::map(gr_l, .remove_substitutions_gr) %>%
      GRangesList()
    return(grl)
  } else if (inherits(grl, "GRanges")) {
    gr <- .remove_substitutions_gr(grl)
    return(gr)
  } else {
    .not_gr_or_grl(grl)
  }
}

#' Removes SNV/MNV variants.
#'
#' @details
#' This function removes SNV/MNV variants for a single GRanges object.
#'
#' @param gr GRanges object
#' @noRd
#' @return A filtered version of the input GRanges object.
#'
.remove_substitutions_gr <- function(gr) {
  snv_f <- .find_substitution(gr)
  gr <- gr[!snv_f]
  return(gr)
}

#' Removes Indel variants.
#'
#' @details
#' This function removes Indel variants for a GRanges/GrangesList object.
#'
#' @param grl GRanges/GRangesList object
#' @noRd
#' @return A filtered version of the input GRanges/GrangesList object.
.remove_indels <- function(grl) {
  if (inherits(grl, "CompressedGRangesList")) {
    gr_l <- as.list(grl)
    grl <- purrr::map(gr_l, .remove_indels_gr) %>%
      GRangesList()
    return(grl)
  } else if (inherits(grl, "GRanges")) {
    gr <- .remove_indels_gr(grl)
    return(gr)
  } else {
    .not_gr_or_grl(grl)
  }
}

#' Removes Indel variants.
#'
#' @details
#' This function removes Indel variants for a single GRanges object.
#'
#' @param gr GRanges object
#' @noRd
#' @return A filtered version of the input GRanges object.
#'
.remove_indels_gr <- function(gr) {
  snv_f <- .find_substitution(gr)
  gr <- gr[snv_f]
  return(gr)
}

#' Identifies SNVs/MNVs
#'
#' @details
#' This function finds SNVs/MNVs for a single GRanges object.
#'
#'
#' @param gr GRanges object
#' @noRd
#' @return A boolean vector. It's TRUE for SNVs/MNVs and FALSE for Indels
#'
.find_substitution <- function(gr) {
  .check_no_multi_alts(gr)
  snv_f <- width(.get_ref(gr)) == width(unlist(.get_alt(gr)))
  return(snv_f)
}

#' Throw an error for when an argument should be a GRanges/GRangesList object.
#'
#' @details
#' This function is called by other functions, when an argument should be a GRanges/GrangesList object, but isn't.
#' It throws an error showing the actual class of the argument.
#'
#' @param arg Argument. Any object that should have been a GRanges/GrangesList, but isn't.
#' @noRd
#' @return An error
#'
.not_gr_or_grl <- function(arg) {
  arg_name <- deparse(substitute(arg))
  arg_class <- class(arg)[[1]]
  stop(stringr::str_c(arg_name, " should be a CompressedGRangesList or GRanges object, instead it is a ", arg_class, " object.
               Please provide an object of the correct class."))
}

#' Throw an error for when an argument should be a GRanges object.
#'
#' @details
#' This function is called by other functions, when an argument should be a GRanges object, but isn't.
#' It throws an error showing the actual class of the argument.
#'
#' @param arg Argument. Any object that should have been a GRanges, but isn't.
#' @noRd
#' @return An error
#'
.not_gr <- function(arg) {
  arg_name <- deparse(substitute(arg))
  arg_class <- class(arg)[[1]]
  stop(stringr::str_c(arg_name, " should be a GRanges object, instead it is a ", arg_class, " object.
               Please provide an object of the correct class."))
}





#' Check that the seqnames of a GRangesList or GRanges object are present in a ref_genome
#'
#' @details
#' This function tests that the variants in a GRangesList or GRanges object have seqnames,
#' that match the supplied ref_genome.
#' If this is not the case an error is thrown, with suggestions how to fix this.
#' This function also check whether variants in the input overlap with the reference genome.
#' If this is not the case, then the wrong reference object might be used.
#'
#' @param grl GRangesList or GRanges object
#' @param ref_genome BSGenome reference genome object
#'
#' @noRd
#'
#' @return Invisibly returns the input grl
#'
.check_chroms <- function(grl, ref_genome) {
  if (inherits(grl, "CompressedGRangesList")) {
    gr <- BiocGenerics::unlist(grl)
  } else if (inherits(grl, "GRanges")) {
    gr <- grl
  } else {
    .not_gr_or_grl(grl)
  }

  # Get seq names
  gr_seqnames <- as.vector(seqnames(gr))
  ref <- BSgenome::getBSgenome(ref_genome)
  ref_seqnames <- seqnames(ref)

  # Check if the gr and reference have the same genome name.
  genome_name_gr <- unique(GenomeInfoDb::genome(gr))
  genome_name_ref <- unique(GenomeInfoDb::genome(ref))
  if (is.na(genome_name_gr) | genome_name_gr != genome_name_ref) {
    stop(paste0(
      "The input GRanges (your vcf data) and the ref_genome do not have the same genome name.\n",
      "This problem is known to occur when you use an outside function to read in vcfs.\n",
      "The input GRanges has the genome name: '", genome_name_gr, "'\n",
      "The ref_genome has the genome name: '", genome_name_ref, "'\n",
      "With `GenomeInfoDb::genome(your input GRanges) = '", genome_name_ref, "'`\n",
      "you can view and change the genome name of your data to that of the ref_genome."
    ), call. = FALSE)
  }

  # Check if there is any overlap in chromosome names
  shared_chroms <- intersect(gr_seqnames, ref_seqnames)
  if (!length(shared_chroms)) {
    stop("The input GRanges and the ref_genome share no seqnames (chromosome names). 
             Do they use the same seqlevelsStyle? An example of how to fix this is show below:
             You can change the seqlevelStyle with: `seqlevelsStyle(grl) = 'UCSC'", call. = FALSE)
  }

  # Check if there seqlevels in the input granges that are not in the reference.
  gr_seqlevels <- levels(seqnames(gr))
  gr_seqlevels_notref <- unique(gr_seqlevels[!gr_seqlevels %in% ref_seqnames])
  if (length(gr_seqlevels_notref)) {
    gr_seqlevels_notref <- stringr::str_c(gr_seqlevels_notref, collapse = ", ")
    stop(stringr::str_c(
      "The following seqlevels (chromosome names) occur in the input GRanges,",
      "but are not present in the ref_genome: ", gr_seqlevels_notref, ". An example of how to fix this is show below.
        First select the chromosomes you want to keep with: 
        `chromosomes = paste0('chr', c(1:22,'X'))`
        You can then remove variants in other chromosomes with: 
        `seqlevels(grl, pruning.mode = 'tidy') = chromosomes`"
    ), call. = FALSE)
  }

  # Check if there are variants in the input granges that don't overlap with the reference
  ref_gr <- GenomicRanges::GRanges(as.vector(seqnames(ref)), IRanges::IRanges(start = 1, end = GenomeInfoDb::seqlengths(ref)))
  hits <- GenomicRanges::findOverlaps(gr, ref_gr)
  if (length(hits)) {
    gr_nomatch <- gr[-S4Vectors::queryHits(hits)]
  } else {
    gr_nomatch <- gr
  }
  nr_nomatch <- length(gr_nomatch)
  if (nr_nomatch) {
    stop(stringr::str_c(
      "There are ", nr_nomatch, " variants that don't overlap with ",
      "the chromosome lengths of the chosen reference genome.\n",
      "Did you select the correct reference genome?"
    ), call. = FALSE)
  }

  invisible(grl)
}


#' Check whether input is a scalar na.
#'
#' @param x Input. Function checks if it is a scalar NA
#'
#' @return Boolean
#' @noRd
#'
.is_na <- function(x) {
  purrr::is_scalar_vector(x) && is.na(x)
}

#' Calculate the amount of lesion segregation for a GRangesList or GRanges object.
#'
#' This function calculates lesion segregation for a GRangesList or GRanges object.
#' Lesion segregation is a large scale Watson versus Crick strand asymmetry caused by
#' many DNA lesions occurring during a single cell cycle.
#' It was first described in Aitken et al., 2020, Nature.
#' See their paper for a more in-depth discussion of this phenomenon.
#' This function can perform three different types of test to calculate lesion segregation.
#' The first method is unique to this package, while the other two were also used by
#' Aitken et al., 2020.
#' The 'binomial' test is based on how often consecutive mutations are on different strands.
#' The 'wald-wolfowitz' test checks if the strands are randomly distributed. It's not known
#' which method is superior.
#' The 'rl20' test looks at run sizes (The number of consecutive mutations on the same strand).
#' This is less susceptible to local strand asymetries and kataegis, but doesn't generate a p-value.
#'
#' The amount of lesion segregation is calculated per GRanges object.
#' The results are then combined in a table.
#'
#' It's possible to calculate the lesion segregation separately per 96 substitution context,
#' when using the binomial test. The results are then automatically added back up together.
#' This can increase sensitivity when a mutational process causes multiple types of base substitutions,
#' which aren’t considered to be on the same strand.
#'
#' When using the rl20 test, this function first calculates the strand runs per chromosome
#' and combines them. It then calculates the smallest set of runs, which together encompass
#' at least 20 percent of the mutations. (This set thus contains the largest runs).
#' The size of the smallest run in this set is the rl20. The genomic span of
#' the runs in this set is also calculated.
#'
#' @param vcf_list GRangesList or GRanges object
#' @param sample_names The name of the sample
#' @param test The statistical test that should be used. Possible values:
#'              * 'binomial' Binomial test based on the number of strand switches. (Default);
#'              * 'wald-wolfowitz' Statistical test that checks if the strands are randomly distributed.;
#'              * 'rl20' Calculates rl20 value and the genomic span of the associated runs set.;
#' @param split_by_type Boolean describing whether the lesion
#' segregation should be calculated for all SNVs together or per 96 substitution context. (Default: FALSE)
#' @param ref_genome BSgenome reference genome object.
#'               Only needed when split_by_type is TRUE with the binomial test
#'               or when using the rl20 test.
#' @param chromosomes The chromosomes that are used. Only needed when using the rl20 test.
#'
#' @return A tibble containing the amount of lesions segregation per sample
#' @importFrom magrittr %>%
#' @seealso
#' \code{\link{plot_lesion_segregation}}
#' @family Lesion_segregation
#' @export
#' @examples
#'
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## To reduce the runtime we take only the first two samples
#' grl <- grl[1:2]
#' ## Set the sample names
#' sample_names <- c("colon1", "colon2")
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Calculate lesion segregation
#' lesion_segretation <- calculate_lesion_segregation(grl, sample_names)
#'
#' ## Calculate lesion segregation per 96 base type
#' lesion_segretation_by_type <- calculate_lesion_segregation(grl, sample_names,
#'   split_by_type = TRUE, ref_genome = ref_genome
#' )
#'
#' ## Calculate lesion segregation using the wald-wolfowitz test.
#' lesion_segregation_wald <- calculate_lesion_segregation(grl,
#'   sample_names,
#'   test = "wald-wolfowitz"
#' )
#'
#' ## Calculate lesion segregation using the rl20.
#' chromosomes <- paste0("chr", c(1:22, "X"))
#' lesion_segregation_rl20 <- calculate_lesion_segregation(grl,
#'   sample_names,
#'   test = "rl20",
#'   ref_genome = ref_genome,
#'   chromosomes = chromosomes
#' )
calculate_lesion_segregation <- function(vcf_list,
                                         sample_names,
                                         test = c("binomial", "wald-wolfowitz", "rl20"),
                                         split_by_type = FALSE,
                                         ref_genome = NA,
                                         chromosomes = NA) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  p.value <- . <- sample_name <- NULL

  # Validate arguments
  test <- match.arg(test)
  if (test != "binomial" & split_by_type) {
    stop("The 'split_by_type' argument can only be used with the binomial test",
      call. = FALSE
    )
  }
  if (length(vcf_list) != length(sample_names)) {
    stop("The vcf_list and the sample_names should be equally long.", call. = FALSE)
  }

  if (.is_na(ref_genome)) {
    if (split_by_type) {
      stop("The ref_genome needs to be set when split_by_type = TRUE", call. = FALSE)
    }
    if (test == "rl20") {
      stop("The ref_genome needs to be set when test == rl20", call. = FALSE)
    }
  }

  if (.is_na(chromosomes) & test == "rl20") {
    stop("The chromosomes need to be set when using test == rl20", call. = FALSE)
  }

  # Turn grl into list.
  if (inherits(vcf_list, "CompressedGRangesList")) {
    vcf_list <- as.list(vcf_list)
  }

  # Perform lesion segregation on each GR
  if (inherits(vcf_list, "list")) {
    strand_tb <- purrr::map2(vcf_list, sample_names, function(gr, sample_name) {
      .calculate_lesion_segregation_gr(gr, sample_name, test, split_by_type, ref_genome, chromosomes)
    }) %>%
      do.call(rbind, .)
  } else if (inherits(vcf_list, "GRanges")) {
    strand_tb <- .calculate_lesion_segregation_gr(
      vcf_list,
      sample_names,
      test,
      split_by_type,
      ref_genome,
      chromosomes
    )
  } else {
    .not_gr_or_grl(vcf_list)
  }

  # Multiple testing correction
  if (test != "rl20") {
    strand_tb <- dplyr::mutate(strand_tb, fdr = p.adjust(p.value, method = "fdr"))
  }

  # Add final columns to output
  strand_tb <- strand_tb %>%
    dplyr::mutate(sample_name = sample_names) %>%
    dplyr::select(sample_name, dplyr::everything())
  return(strand_tb)
}

#' Calculate the amount of lesion segregation for a singe GRanges object.
#'
#' @param gr GRanges object
#' @param sample_name The name of the sample
#' @param test The statistical test that should be used. Possible values:
#'              * 'binomial' Binomial test based on the number of strand switches. (Default);
#'              * 'wald-wolfowitz' Statistical test that checks if the strands are randomly
#'              distributed.;
#'              * 'rl20' Calculates rl20 value and the genomic span of the associated runs set.;
#' @param split_by_type Boolean describing whether the lesion
#' segregation should be calculated for all SNVs together or per 96 substitution context.
#' @param ref_genome BSgenome reference genome object.
#'               Only needed when split_by_type is TRUE with the binomial test
#'               or when using the rl20 test.
#' @param chromosomes The chromosomes that are used. Only needed when using the rl20 test.
#'
#' @return A tibble containing the amount of lesions segregation for a single sample
#' @noRd
#'
.calculate_lesion_segregation_gr <- function(gr,
                                             sample_name = "sample",
                                             test = c("binomial", "wald-wolfowitz", "rl20"),
                                             split_by_type = FALSE,
                                             ref_genome = NA,
                                             chromosomes = NA) {


  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  genome_span <- genome_size <- NULL

  # Check if mutations are present.
  if (!length(gr)) {
    message(paste0(
      "No mutations present in sample: ", sample_name,
      "\n Returning NA"
    ))
    return(NA)
  }

  # Get strand info
  gr <- .get_strandedness_gr(gr)
  tb <- .get_strandedness_tb(gr)

  if (test == "binomial") {
    # Perform analysis per base substitution type
    if (split_by_type) {


      # Split gr according to the 96 substitution context.
      cnd <- tryCatch(suppressWarnings({
        GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
      }),
      error = function(cnd) cnd
      )
      if (inherits(cnd, "error")) {
        message(paste0(
          "Could not change seqlevelstyle in sample: ", sample_name, ".",
          "\n Returning NA"
        ))
        return(NA)
      }
      .check_chroms(gr, ref_genome)
      type_context <- type_context(gr, ref_genome)
      full_context <- paste0(
        substr(type_context$context, 1, 1),
        "[", type_context$types, "]",
        substr(type_context$context, 3, 3)
      )
      tb_l <- split(tb, full_context)

      # Calculate strand switches for each of the 96 substitutions.
      res_l <- purrr::map(tb_l, .calculate_strand_switches)
      x <- purrr::map(res_l, "x") %>%
        unlist() %>%
        sum()
      n <- purrr::map(res_l, "n") %>%
        unlist() %>%
        sum()
      res <- list("x" = x, "n" = n)
    } else {
      # Calculate strand switches
      res <- .calculate_strand_switches(tb)
    }

    # Check if mutations are present
    if (res$n == 0) {
      message(paste0(
        "No multiple mutations in one chromosome with context present in sample: ", sample_name,
        "\n Returning NA"
      ))
      return(NA)
    }

    # Calculate if the number of strand switches is significantly different from expected.
    res <- binom.test(x = res$x, n = res$n, p = 0.5)

    # Add all results together in a tibble
    stat_tb <- tibble::tibble(
      p.value = res$p.value,
      fraction_strand_switches = res$estimate,
      conf_low = res$conf.int[[1]],
      conf_high = res$conf.int[[2]],
      nr_strand_switches = res$statistic,
      max_possible_switches = res$parameter
    )
  } else if (test == "wald-wolfowitz") {
    # calculate if there is a significant deviation using the wald_wolfowitz_test
    wolfowitz <- .wald_wolfowitz_test(tb$strand)
    stat_tb <- tibble::tibble(
      p.value = wolfowitz$p,
      sd = wolfowitz$sd,
      nr_total_runs = wolfowitz$runs_total
    )
  } else if (test == "rl20") {

    # Calculate rl20 and genomic span
    res <- .rl20_gspan(tb)

    # Add total size of genome to calculate fraction.
    ref_genome <- BSgenome::getBSgenome(ref_genome)
    stat_tb <- res %>%
      dplyr::mutate(
        genome_size = sum(GenomeInfoDb::seqlengths(ref_genome)[chromosomes]),
        fraction_span = genome_span / genome_size
      )
  }

  return(stat_tb)
}

#' Determine the strands of a GRanges object
#'
#' @param gr A GRanges object
#'
#' @return A GRanges object where the strands have been set.
#' @noRd
#'
.get_strandedness_gr <- function(gr) {
  .check_no_indels(gr)
  strand(gr) <- ifelse(as.vector(.get_ref(gr)) %in% c("C", "T"), "+", "-")

  if (length(gr)) {
    GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::seqlevelsInUse(gr)
  }
  return(gr)
}

#' Convert a GRanges object with strand info to a tibble
#'
#' @param gr A GRanges object where the strands have been set.
#'
#' @return A tibble with strand information
#' @importFrom magrittr %>%
#' @noRd
#'
.get_strandedness_tb <- function(gr) {
  tb <- as.data.frame(gr) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      strand = droplevels(strand),
      y = dplyr::recode(strand, "+" = 1, "-" = 0),
      start_mb = start / 1000000
    )
  return(tb)
}

#' Calculate the total number of variants and strand switches.
#'
#' @param tb A tibble with strand information
#'
#' @return A list containing the total number of variants and the number of strand switches
#' @importFrom magrittr %>%
#' @noRd
#'
.calculate_strand_switches <- function(tb) {
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  . <- NULL

  strands_l <- split(tb$strand, tb$seqnames)
  switches <- purrr::map(strands_l, .calculate_strand_switch) %>%
    do.call("c", .)
  res <- list("x" = sum(switches), "n" = length(switches))
  return(res)
}

#' Calculate the number of strand switches in a single chromosome.
#'
#' @param strand A vector containing the strand of variants
#'
#' @return A boolean vector describing for each variant
#' if it switched strands with the previous variant.
#' @noRd
#' @importFrom magrittr %>%
#'
.calculate_strand_switch <- function(strand) {
  switches <- strand != dplyr::lead(strand)
  switches <- switches %>%
    stats::na.omit() %>%
    as.vector()
  return(switches)
}


#' Perform the wald_wofowitz test for strands.
#'
#' This statistical test, tests whether each element in the sequence is
#' independently drawn from the same distribution.
#'
#' @param strands A vector of strands
#'
#' @return a p value
#'
#' @noRd
#'
.wald_wolfowitz_test <- function(strands) {

  # Remove factor
  strands <- as.character(strands)

  # Determine sizes
  n1 <- sum(strands == "+")
  n2 <- sum(strands == "-")
  n <- n1 + n2

  # Determine number of + and - runs
  runs <- rle(strands)
  r1 <- length(runs$lengths[runs$values == "+"])
  r2 <- length(runs$lengths[runs$values == "-"])

  # Calculate total number of runs
  runs_total <- r1 + r2

  # Calculate mean
  mean_val <- 2 * n1 * n2 / (n) + 1

  # Calculate variance and sd
  variance <- (mean_val - 1) * (mean_val - 2) / (n - 1)
  sd <- sqrt(variance)

  # Calculate p value
  p <- stats::pnorm((runs_total - mean_val) / sd)

  # Make two-sided
  p <- 2 * min(p, 1 - p)

  return(list("p" = p, "sd" = sd, "runs_total" = runs_total))
}

#' Calculate rl20 and genomic span
#'
#' This function calculates the strand runs per chromosome
#' and combines them. It then calculates the rl20. This is done by
#' calculating the smallest set of runs, which together encompass
#' at least 20% of the mutations. (This set thus contains the largest runs).
#' The size of the smallest run in this set is the rl20. The genomic span of
#' the runs in this set is also calculated.
#'
#' @param tb A tibble with strand information
#'
#' @return A list containing the rl20 and the genomic span
#' @noRd
#'
.rl20_gspan <- function(tb) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  . <- NULL

  # Number of variants
  n <- nrow(tb)

  # Remove factor
  tb <- dplyr::mutate(tb, strand = as.character(strand))

  # Determine runs per chromosome
  tb_l <- split(tb, tb$seqnames)
  runs_l <- purrr::map(tb_l, ~ rle(.x$strand))

  # Combine run lengths
  run_lengths <- purrr::map(runs_l, "lengths") %>%
    do.call(c, .) %>%
    magrittr::set_names(NULL)


  # Sort runs on size
  order_i <- order(run_lengths, decreasing = TRUE)
  sort_lengths <- run_lengths[order_i]

  # Determine set of runs that together encompass 20% of mutations.
  not_big_enough_set_size <- sum(cumsum(sort_lengths) < 0.2 * n) # Just less than 20%
  set_size <- not_big_enough_set_size + 1 #+1 to get over 20%
  set <- sort_lengths[seq_len(set_size)]

  # Get shortest run from 20% set
  rl20 <- set[set_size]

  # The next section will determine the genomic span.
  # To do this we need to get the position
  # of variants back from the runs.

  # Determine original index of the runs in set.
  order_i_set <- order_i[seq_len(set_size)]

  # Get the mutation indices from the runs
  cumsum_runs <- cumsum(run_lengths)

  # Determine index of the first mutation of the runs in set.
  # The index of the last mutation in the previous run is used for this.
  # If the first run is in the rl20, then it has to be set separately.
  first_run_i <- which(order_i_set == 1)
  if (length(first_run_i)){
    order_i_set[first_run_i] <- 2
    start_i <- cumsum_runs[order_i_set - 1] + 1
    start_i[first_run_i] <- 1
  } else{
    start_i <- cumsum_runs[order_i_set - 1] + 1
  }
  # Determine index of the last mutation of the runs in set.
  end_i <- cumsum_runs[order_i_set]

  # Use the indices to get the genomic positions and calculate
  # the genomic span.
  genomic_spans <- tb$end[end_i] - tb$start[start_i]
  genomic_span <- sum(genomic_spans)

  res <- tibble::tibble("rl20" = rl20, "genome_span" = genomic_span)
  return(res)
}

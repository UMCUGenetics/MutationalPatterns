#' Plot the strands of variants to show lesion segregation
#'
#' The strands of variants in a GRanges object is plotted.
#' This way the presence of any lesion segregation is visualized.
#' The function can plot either a single or multiple samples.
#' Per chromosome, the ratio of the mutations on the chromosomal strands is
#' visualised by a line. The position of this line is calculated as the mean of
#' the "+" and "-" strand, where "+" equals 1 and "-" equals 0. In other words:
#' this line lies between the two strands if the mutations are equally
#' distributed between them, and approaches a strand if the majority of
#' mutations on a chromosome lie on that strand.
#' 
#' @param vcf GRanges or RGrangesList object.
#' @param per_chrom Boolean. Determines whether to create a separate plot per chromosome.
#' @param sample_name Name of the sample. Is used as the title of the plot.
#' Not very useful if you have more than one sample.
#' @param min_muts_mean Integer. The minimum of mutations, required for the mean strand 
#' of a chromosome to be calculated.
#' @param chromosomes Character vector. Determines chromosomes to be used and their order.
#' @param subsample Double between 0 and 1. Subsamples the amount of mutations to create
#' a plot with less dots. Such a plot is easier to modify in a vector program like illustrator.
#' (default: NA)
#'
#' @return ggplot2 object
#' @export
#' @seealso
#' \code{\link{calculate_lesion_segregation}}
#' @family Lesion_segregation
#' @examples
#'
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#' 
#' ## Plot lesion segregation
#' plot_lesion_segregation(grl[1:3])
#' 
#' ## Select a single GRanges object to plot.
#' gr <- grl[[1]]
#'
#' ## Plot lesion segregation for a single sample. 
#' ## Also add a title to the plot.
#' plot_lesion_segregation(gr, sample_name = "Colon1")
#'
#' ## Plot lesion segregation per chromosome.
#' ## We here store the results in a list.
#' figure_l = plot_lesion_segregation(gr, per_chrom = TRUE, sample_name = "Colon1")
#' 
#' ## Plot specific chromosomes in a user specified order
#' plot_lesion_segregation(grl[1:3], chromosomes = c(2,3))
#' 
#' ## Subsample the mutations, so less points are plotted.
#' plot_lesion_segregation(grl[1:3], subsample = 0.2)
#' 
plot_lesion_segregation <- function(vcf, 
                                    per_chrom = FALSE, 
                                    sample_name = NA, 
                                    min_muts_mean = 10,
                                    chromosomes = NA,
                                    subsample = NA) {

  # Argument name is vcf, because the function previously only worked on a single GRanges.
  # Because the function now also works on a GRangesList, the name vcf_list is used internally.
  vcf_list <- vcf
  
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  max_pos <- start_mb <- notused <- n <- y <- NULL

  # Convert list to grl if necessary
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  }
  # Genome is set to NULL to ensure seqlevels can be changed.
  GenomeInfoDb::genome(vcf_list) <- NA
  GenomeInfoDb::seqlevelsStyle(vcf_list) <- "NCBI" # This takes less space when plotting
  
  
  # Subsample the amount of mutations if necessary
  if (!is.na(subsample)){
    if (inherits(vcf_list, "CompressedGRangesList")){
      vcf_list <- purrr::map(as.list(vcf_list), .subsample_granges, subsample) %>% 
        GenomicRanges::GRangesList()
    } else if (inherits(vcf_list, "GRanges")) {
      vcf_list <- .subsample_granges(vcf_list, subsample)
    } else {
      .not_gr_or_grl(vcf_list)
    }
  }

  
  # get strandedness
  if (inherits(vcf_list, "CompressedGRangesList")){
    tb <- purrr::map(as.list(vcf_list), function(vcf) {
      gr <- .get_strandedness_gr(vcf)
      tb <- .get_strandedness_tb(gr)
      return(tb)
    }) %>% 
      dplyr::bind_rows(.id = "sample")
    sample_names <- names(vcf_list)
    nr_samples <- length(vcf_list)
    max_nr_mutations <- max(S4Vectors::elementNROWS(vcf_list))
  } else if (inherits(vcf_list, "GRanges")) {
    gr <- .get_strandedness_gr(vcf_list)
    tb <- .get_strandedness_tb(gr)
    tb$sample <- "MySample"
    sample_names <- "MySample"
    nr_samples <- 1
    max_nr_mutations <- length(vcf_list)
  } else {
    .not_gr_or_grl(vcf_list)
  }
  
  # Turn the sample column into a factor.
  # This ensures the correct order of the samples in the plot.
  tb <- dplyr::mutate(tb, sample = factor(sample, levels = unique(sample)))
  
  
  ## set chrom order
  if (!.is_na(chromosomes)) {
    chromosomes <- as.character(chromosomes)
    withCallingHandlers(
      {
        GenomeInfoDb::seqlevelsStyle(chromosomes) = "NCBI"
      },
      warning = function(w) {
        if (grepl("found more than one best sequence renaming map compatible with seqname style", 
                  conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
    
    tb <- tb %>%
      dplyr::filter(seqnames %in% chromosomes) %>%
      dplyr::mutate(seqnames = factor(seqnames, levels = chromosomes))
  } else{
    chromosomes <- GenomeInfoDb::seqlevelsInUse(vcf_list)
    if (!length(chromosomes)){
      chromosomes <- GenomeInfoDb::seqlevels(vcf_list)
    }
    tb$seqnames <- factor(tb$seqnames, levels = chromosomes)
  }
  
  # Create limit to ensure that the entire chromosomes are plotted, 
  # even when mutations don't span the entire chromosome.
  chroms <- levels(tb$seqnames)
  nr_chroms <- length(chroms)
  chrom_lengths <- GenomeInfoDb::seqlengths(vcf_list)[chroms]
  
  chr_starts <- tibble::tibble("start_mb" = rep(1, nr_chroms),
                              "seqnames" = chroms)
  chr_ends <- tibble::tibble("start_mb" = chrom_lengths/1000000,
                             'seqnames' = chroms)
  chr_limits <- rbind(chr_starts, chr_ends)
  
  # Per sample limits
  tb_limits <- tibble::tibble("sample" = factor(rep(sample_names, each = 2 * nr_chroms), levels = levels(tb$sample)),
                          "start_mb" = rep(chr_limits$start_mb, nr_samples),
                          "seqnames" = factor(rep(chr_limits$seqnames, nr_samples), levels = chroms),
                          "y" = 0.5)

  
  # calculate mean strand per chomosome per sample
  tb_mean <- tb %>% 
    dplyr::group_by(seqnames, sample) %>% 
    dplyr::summarize(y = mean(y), n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::ungroup() %>% 
    dplyr::right_join(tb_limits[,c("sample", "seqnames", "start_mb")], 
                                by = c("seqnames", "sample")) %>% 
    dplyr::filter(n >= min_muts_mean) %>% 
    dplyr::mutate(sample = factor(sample, levels = levels(tb$sample)))

  # Set point_sizes
  point_size <- 200 / max_nr_mutations
  if (per_chrom == TRUE) {
    point_size <- point_size * 5
  }
  if (point_size > 2) {
    point_size <- 2
  } else if (point_size < 0.02) {
    point_size <- 0.02
  }

  # Create plots
  if (per_chrom == FALSE) {
    fig <- .plot_lesion_segregation_gg(tb, tb_mean, tb_limits, point_size, sample_name)
    return(fig)
  } else {
    tb_l <- split(tb, tb$seqnames)
    tb_mean_l <- split(tb_mean, tb_mean$seqnames)
    tb_limits_l <- split(tb_limits, tb_limits$seqnames)
    fig_l <- mapply(.plot_lesion_segregation_gg,
      tb_l, tb_mean_l, tb_limits_l,
      MoreArgs = list("point_size" = point_size, "sample_name" = sample_name),
      SIMPLIFY = FALSE
    )
    return(fig_l)
  }
}


#' Subsample granges to a smaller number of mutations
#'
#' @param vcf GRanges object
#' @param subsample Double between 0 and 1. Subsamples the amount of mutations to create
#' a plot with less dots. Such a plot is easier to modify in a vector program like illustrator.
#' 
#' @return GRanges object
#' @noRd
#'
.subsample_granges = function(vcf, subsample){
  nr_muts_kept <- ceiling(subsample * length(vcf))
  muts_kept_i <- sample.int(length(vcf), nr_muts_kept, replace = FALSE)
  vcf <- vcf[muts_kept_i,]
  return(vcf)
}

#' Plot the strands of variants to show lesion segregation
#'
#' This is a helper function for 'plot_lesion_segregation'.
#' It performs the actual plotting.
#'
#' @param tb A tibble with strand information of variants
#' @param tb_mean A tibble with the mean strand per chromosome
#' per sample.
#' @param tb_limits tibble describing the chromosome limits.
#' This ensures entire chromosomes are plotted,
#' instead of just the part with variants
#' @param point_size Scalar describing the point size of the plot
#' @param sample_name Name of the sample
#'
#' @return ggplot2 object
#'
#' @import ggplot2
#' @noRd
#'
.plot_lesion_segregation_gg <- function(tb, tb_mean, tb_limits, point_size, sample_name) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  y <- start_mb <- NULL

  if (.is_na(sample_name)) {
    my_labs <- labs(y = "Strand", x = "Coordinate (mb)")
  } else {
    my_labs <- labs(y = "Strand", x = "Coordinate (mb)", title = sample_name)
  }


  # Plot strandedness
  fig <- ggplot(mapping = aes(y = y, x = start_mb, color = y)) +
    geom_line(data = tb_mean, size = 1.5) +
    geom_jitter(data = tb, width = 0.05, height = 0.1, size = point_size) +
    facet_grid(sample ~ seqnames, scales = "free_x", space = "free_x") +
    geom_blank(data = tb_limits) +
    scale_y_continuous(breaks = c(0, 1), labels = c("-", "+"), limits = c(-1, 2)) +
    my_labs +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 16),
          panel.grid = element_blank(),
          panel.grid.major.y = element_line(),
          panel.border = element_blank(),
          legend.position = 'none') +
    scale_color_gradientn(colors = c('#446DF6', 'grey90', '#FFBC63'),limits = c(0,1))
  
  return(fig)
}

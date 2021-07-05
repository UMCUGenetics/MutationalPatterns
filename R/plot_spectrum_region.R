#' Plot point mutation spectrum per genomic region
#'
#' A spectrum similar to the one from 'plot_spectrum()' is plotted.
#' However the spectrum is plotted separately per genomic region.
#' As input it takes a 'type_occurrences' matrix that was calculated per genomic region.
#' To get a 'type_occurrences' matrix per region,
#' first use the 'split_muts_region()' function on a GR or GRangesList.
#' Then use the 'mut_type_occurrences' as you would normally.
#' The by, colors and legend argument work the same as in 'plot_spectrum()'.
#'
#' The y-axis can be plotted with three different modes. With
#' 'relative_sample_feature', the number of variants will be shown divided by
#' the total number of variants in that sample and genomic region. This is
#' generally the most useful, because it allows you to compare the spectra off
#' different regions. When you use 'relative_sample', the number of variants
#' will be shown divided by the total number of variants in that sample. This
#' can be useful when you want to compare the number of mutations between
#' regions. Finally, when you use 'absolute', the absolute mutation numbers are
#' shown. This can be useful when you want to compare the mutation load between
#' different groups of samples.
#'
#'
#' @param type_occurrences Type occurrences matrix
#' @param by Optional grouping variable
#' @param mode  The y-axis plotting mode.
#'              * 'relative_sample', the number of variants will be shown
#'              divided by the total number of variants in that sample;
#'              * 'relative_sample_feature', the number of variants will be shown
#'              divided by the total number of variants in that sample and genomic region (Default);
#'              * 'absolute' The absolute number of mutations is shown;
#' @param indv_points Whether to plot the individual samples
#' as points, default = FALSE
#' @param error_bars The type of error bars to plot.
#'              * '95%_CI' for 95% Confidence intervals (Default);
#'              * 'stdev' for standard deviations;
#'              * 'SEM' for the standard error of the mean (NOT recommended);
#'              * 'none' Do not plot any error bars;
#' @param colors Optional color vector with 7 values
#' @param legend Plot legend, default = TRUE
#' @param condensed More condensed plotting format. Default = F.
#' @return Spectrum plot by genomic region
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' ## See the 'split_muts_region()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/grl_split_region.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#'
#' ## Get the type occurrences for all VCF objects.
#' type_occurrences <- mut_type_occurrences(grl, ref_genome)
#'
#' ## Plot the relative point mutation spectrum per genomic region
#' plot_spectrum_region(type_occurrences)
#'
#' ## Include the individual sample points
#' plot_spectrum_region(type_occurrences, indv_points = TRUE)
#'
#' ## Plot the relative point mutation spectrum per genomic region,
#' ## but normalize only for the samples
#' plot_spectrum_region(type_occurrences, mode = "relative_sample")
#'
#' ## Plot the absolute point mutation spectrum per genomic region
#' plot_spectrum_region(type_occurrences, mode = "absolute")
#'
#' ## Plot the point mutations spectrum with different error bars
#' plot_spectrum_region(type_occurrences, error_bars = "stdev")
#'
#' ## Plot the relative point mutation spectrum per sample type and per genomic region
#' ## Determine tissue names
#' tissue <- c(
#'   "colon", "colon", "colon",
#'   "intestine", "intestine", "intestine",
#'   "liver", "liver", "liver"
#' )
#' plot_spectrum_region(type_occurrences, by = tissue)
#'
#' ## Plot the relative point mutation spectrum per individual sample and per genomic region
#' ## Determine sample names
#' sample_names <- c(
#'   "colon1", "colon2", "colon3",
#'   "intestine1", "intestine2", "intestine3",
#'   "liver1", "liver2", "liver3"
#' )
#'
#' plot_spectrum_region(type_occurrences, by = sample_names, error_bars = "none")
#' 
#' ## Plot it in a more condensed manner, 
#' ## which is is ideal for publications.
#' plot_spectrum_region(type_occurrences, 
#' by = sample_names, 
#' error_bars = "none", 
#' condensed = TRUE)
#' 
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_type_occurrences}},
#' \code{\link{plot_spectrum}},
#' \code{\link{split_muts_region}}
#' @family genomic_regions
#'
#'
plot_spectrum_region <- function(type_occurrences,
                                 by = NA,
                                 mode = c("relative_sample_feature", "relative_sample", "absolute"),
                                 indv_points = FALSE,
                                 error_bars = c("95%_CI", "stdev", "SEM", "none"),
                                 colors = NULL,
                                 legend = TRUE,
                                 condensed = FALSE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  `C>T at CpG` <- `C>T other` <- type <- amount <- stdev <- tot_muts <- lower <- upper <- freq <- NULL
  total_indv <- sem <- error_95 <- NULL

  # Check arguments
  if (is.null(colors)) {
    colors <- COLORS6
  }
  mode <- match.arg(mode)
  error_bars <- match.arg(error_bars)

  row_names <- rownames(type_occurrences)
  max_dots_in_name <- row_names %>%
    stringr::str_count("\\.") %>%
    max()
  if (max_dots_in_name > 1) {
    stop("The row names of the type_occurrences dataframe too many dots. 
             There should only be a dot in between the sample name and the type")
  }

  # Get sample names and features
  sample_names <- stringr::str_remove(row_names, "\\..*")
  feature <- stringr::str_remove(row_names, ".*\\.")
  feature <- factor(feature, levels = unique(feature))

  # Remove CpG split
  type_occurrences <- type_occurrences %>%
    dplyr::select(-`C>T at CpG`, -`C>T other`)


  # Create long tb
  tb_per_sample <- type_occurrences %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sample = sample_names, feature = feature) %>%
    tidyr::gather(key = "type", value = "amount", -sample, -feature)

  # Normalize depending on mode
  if (mode == "relative_sample") {
    tb_per_sample <- tb_per_sample %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(freq = amount / sum(amount)) %>%
      dplyr::ungroup()
    y_lab <- "Relative contribution"
  } else if (mode == "relative_sample_feature") {
    tb_per_sample <- tb_per_sample %>%
      dplyr::group_by(sample, feature) %>%
      dplyr::mutate(freq = amount / sum(amount)) %>%
      dplyr::ungroup()
    y_lab <- "Relative contribution"
  } else if (mode == "absolute") {
    tb_per_sample <- tb_per_sample %>%
      dplyr::mutate(freq = amount)
    y_lab <- "Contribution"
  }
  tb_per_sample <- dplyr::mutate(tb_per_sample, freq = ifelse(is.nan(freq), 0, freq))

  # Add sample grouping
  if (.is_na(by)) {
    by <- "all"
  }
  tb_by <- tibble::tibble(
    "sample" = unique(tb_per_sample$sample),
    "by" = by
  )
  tb_per_sample <- tb_per_sample %>%
    dplyr::left_join(tb_by, by = "sample")

  # Combine samples based on sample grouping
  tb <- tb_per_sample %>%
    dplyr::mutate(by = factor(by, levels = unique(by))) %>% 
    dplyr::group_by(by, feature, type) %>%
    dplyr::summarise(
      stdev = sd(freq),
      freq = mean(freq),
      amount = sum(amount),
      total_indv = dplyr::n(),
      .groups = "drop_last"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sem = stdev / sqrt(total_indv),
      error_95 = ifelse(total_indv > 1, qt(0.975, df = total_indv - 1) * sem, NA)
    )

  # Count nr muts per sample group
  tot_muts_tb <- tb %>%
    dplyr::group_by(by) %>%
    dplyr::summarise(tot_muts = sum(amount))

  # Create facet labels
  facet_labs_y <- stringr::str_c(tot_muts_tb$by, " (n = ", tot_muts_tb$tot_muts, ")")
  names(facet_labs_y) <- tot_muts_tb$by

  # Change plotting parameters based on whether plot should be condensed.
  if (condensed == TRUE) {
    spacing <- 0
  } else {
    spacing <- 0.5
  }
  
  # Create figure
  # Suppress warning about using alpha.
  withCallingHandlers(
    {
      fig <- ggplot(tb, aes(x = type, 
                            y = freq, 
                            fill = type, 
                            alpha = feature)) +
        geom_bar(stat = "identity", position = "dodge", colour = "black", cex = 0.5) +
        facet_grid(. ~ by, labeller = labeller(by = facet_labs_y)) +
        scale_fill_manual(values = colors) +
        scale_alpha_discrete(range = c(1, 0.4)) +
        labs(y = y_lab, x = "") +
        scale_x_discrete(breaks = NULL) +
        theme_bw() +
        theme(panel.spacing.x = unit(spacing, "lines"))
    },
    warning = function(w) {
      if (grepl("Using alpha for a discrete variable is not advised.", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )

  # Add individual points
  if (indv_points == TRUE) {
    # Add total_mutations column, which is necessary for faceting later
    fig <- fig +
      geom_point(
        data = tb_per_sample, aes(y = freq), colour = "grey23", shape = 21,
        position = position_jitterdodge(dodge.width = 1, jitter.width = 0.3)
      )
  }

  # Add errorbars
  if (sum(is.na(tb$stdev)) != 0 & error_bars != "none") {
    warning("No error bars can be plotted, because there is only one sample per mutation spectrum.
              Use the argument: `error_bars = 'none'`, if you want to avoid this warning.",
      call. = FALSE
    )
  } else {
    if (error_bars == "stdev") {
      fig <- fig + geom_errorbar(aes(
        ymin = freq - stdev,
        ymax = freq + stdev
      ), position = "dodge")
    } else if (error_bars == "95%_CI") {
      fig <- fig + geom_errorbar(aes(
        ymin = freq - error_95,
        ymax = freq + error_95
      ), position = "dodge")
    } else if (error_bars == "SEM") {
      fig <- fig + geom_errorbar(aes(
        ymin = freq - sem,
        ymax = freq + sem
      ), position = "dodge")
    }
  }

  # Remove legend if required
  if (legend == FALSE) {
    fig <- fig + guides(fill = "none", alpha = "none")
  }

  return(fig)
}

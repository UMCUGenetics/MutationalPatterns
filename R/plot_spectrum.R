#' Plot point mutation spectrum
#'
#' @param type_occurrences Type occurrences matrix
#' @param CT Distinction between C>T at CpG and C>T at other
#' sites, default = FALSE
#' @param by Optional grouping variable
#' @param indv_points Whether to plot the individual samples
#' as points, default = FALSE
#' @param error_bars The type of error bars to plot.
#'              * '95%_CI' for 95% Confidence intervals (default);
#'              * 'stdev' for standard deviations;
#'              * 'SEM' for the standard error of the mean (NOT recommended);
#'              * 'none' Do not plot any error bars;
#' @param colors Optional color vector with 7 values
#' @param legend Plot legend, default = TRUE
#' @param condensed More condensed plotting format. Default = F.
#' @return Spectrum plot
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Get the type occurrences for all VCF objects.
#' type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
#'
#' ## Plot the point mutation spectrum over all samples
#' plot_spectrum(type_occurrences)
#'
#' ## Or with distinction of C>T at CpG sites
#' plot_spectrum(type_occurrences, CT = TRUE)
#'
#' ## You can also include individual sample points.
#' plot_spectrum(type_occurrences, CT = TRUE, indv_points = TRUE)
#'
#' ## You can also change the type of error bars
#' plot_spectrum(type_occurrences, error_bars = "stdev")
#'
#' ## Or plot spectrum per tissue
#' tissue <- c(
#'   "colon", "colon", "colon",
#'   "intestine", "intestine", "intestine",
#'   "liver", "liver", "liver"
#' )
#'
#' plot_spectrum(type_occurrences, by = tissue, CT = TRUE)
#'
#' ## Or plot the spectrum per sample. Error bars are set to 'none', because they can't be plotted.
#' plot_spectrum(type_occurrences, by = names(vcfs), CT = TRUE, error_bars = "none")
#'
#' ## Plot it in a more condensed manner, 
#' ## which is is ideal for publications.
#' plot_spectrum(type_occurrences, 
#' by = names(vcfs), 
#' CT = TRUE, 
#' error_bars = "none",
#' condensed = TRUE)
#'
#' ## You can also set custom colors.
#' my_colors <- c(
#'   "pink", "orange", "blue", "lightblue",
#'   "green", "red", "purple"
#' )
#'
#' ## And use them in a plot.
#' plot_spectrum(type_occurrences,
#'   CT = TRUE,
#'   legend = TRUE,
#'   colors = my_colors
#' )
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_type_occurrences}}
#'
#' @export

plot_spectrum <- function(type_occurrences, 
                          CT = FALSE, 
                          by = NA, 
                          indv_points = FALSE,
                          error_bars = c("95%_CI", "stdev", "SEM", "none"), 
                          colors = NA, 
                          legend = TRUE,
                          condensed = FALSE) {
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  value <- nmuts <- sub_type <- variable <- error_pos <- stdev <- total_mutations <- NULL
  x <- total_individuals <- sem <- error_95 <- NULL

  # Match argument
  error_bars <- match.arg(error_bars)

  # If colors parameter not provided, set to default colors
  if (.is_na(colors)) {
    colors <- COLORS7
  }

  # Check color vector length
  if (length(colors) != 7) {
    stop("Colors parameter: supply color vector with length 7")
  }

  # Distinction between C>T at CpG or not
  if (CT == FALSE) {
    type_occurrences <- type_occurrences[, seq_len(6)]
  } else {
    type_occurrences <- type_occurrences[, c(1, 2, 8, 7, 4, 5, 6)]
  }

  # If grouping variable not provided, set to "all"
  if (.is_na(by)) {
    by <- "all"
  }

  # Reshape the type_occurences for the plotting
  tb_per_sample <- type_occurrences %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::mutate(by = by) %>% # Add user defined grouping
    tidyr::pivot_longer(c(-sample, -by), names_to = "variable", values_to = "nmuts") %>% # Make long format
    dplyr::group_by(sample) %>%
    dplyr::mutate(value = nmuts / sum(nmuts)) %>% # Calculate relative values
    dplyr::ungroup() %>%
    dplyr::mutate(
      sub_type = stringr::str_remove(variable, " .*"),
      variable = factor(variable, levels = unique(variable))
    )
  
  # Summarise per group and mutation type
  tb <- tb_per_sample %>%
    dplyr::mutate(by = factor(by, levels = unique(by))) %>% 
    dplyr::group_by(by, variable) %>%
    dplyr::summarise(
      sub_type = sub_type[[1]], mean = mean(value), stdev = stats::sd(value),
      total_individuals = sum(value), total_mutations = sum(nmuts)
    ) %>%
    dplyr::mutate(total_individuals = sum(total_individuals), total_mutations = sum(total_mutations)) %>%
    dplyr::mutate( # Calc 95% CI and sem
      sem = stdev / sqrt(total_individuals),
      error_95 = ifelse(total_individuals > 1, qt(0.975, df = total_individuals - 1) * sem, NA)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate( # Make pretty and add subtypes
      total_mutations = prettyNum(total_mutations, big.mark = ","),
      total_mutations = paste("No. mutations = ", total_mutations),
      error_pos = mean
    )

  # Change some settings based on whether CT should be plotted separately.
  if (CT == FALSE) {
    # Define colors for plotting
    colors <- colors[c(1, 2, 3, 5, 6, 7)]
  } # C>T stacked bar (distinction between CpG sites and other)
  else {
    # Adjust positioning of error bars for stacked bars
    # mean of C>T at CpG should be plus the mean of C>T at other
    CpG <- which(tb$variable == "C>T at CpG")
    other <- which(tb$variable == "C>T other")
    tb$error_pos[CpG] <- tb$error_pos[other] + tb$error_pos[CpG]

    # Value of the individual sample points also needs to be adjusted.
    CpG <- which(tb_per_sample$variable == "C>T at CpG")
    other <- which(tb_per_sample$variable == "C>T other")
    tb_per_sample$value[CpG] <- tb_per_sample$value[other] + tb_per_sample$value[CpG]
  }

  # Change plotting parameters based on whether plot should be condensed.
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  } else {
    width <- 0.9
    spacing <- 0.5
  }
  
  # Make barplot
  plot <- ggplot(data = tb, aes(
    x = sub_type,
    y = mean,
    fill = variable,
    group = sub_type,
    width = width
  )) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, name = "Point mutation type") +
    theme_bw() +
    xlab("") +
    ylab("Relative contribution") +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = unit(spacing, "lines")
    )

  # Add individual points
  if (indv_points == TRUE) {
    # Add total_mutations column, which is necessary for faceting later
    tb_per_sample <- dplyr::left_join(tb_per_sample,
      tb[, c("by", "variable", "total_mutations")],
      by = c("by", "variable")
    )
    plot <- plot +
      geom_jitter(
        data = tb_per_sample, aes(y = value),
        height = 0, width = 0.3, shape = 21, colour = "grey23"
      )
  }

  # Add error bars
  if (sum(is.na(tb$stdev)) > 0 & error_bars != "none") {
    warning("No error bars can be plotted, because there is only one sample per mutation spectrum.
              Use the argument: `error_bars = 'none'`, if you want to avoid this warning.",
      call. = FALSE
    )
  }
  else {
    if (error_bars == "stdev") {
      plot <- plot + geom_errorbar(aes(
        ymin = error_pos - stdev,
        ymax = error_pos + stdev
      ), width = 0.2)
    } else if (error_bars == "95%_CI") {
      plot <- plot + geom_errorbar(aes(
        ymin = error_pos - error_95,
        ymax = error_pos + error_95
      ), width = 0.2)
    } else if (error_bars == "SEM") {
      plot <- plot + geom_errorbar(aes(
        ymin = error_pos - sem,
        ymax = error_pos + sem
      ), width = 0.2)
    }
  }


  # Perform facetting
  if (length(by) == 1) {
    plot <- plot + facet_wrap(~total_mutations)
  } else {
    plot <- plot + facet_wrap(by ~ total_mutations)
  }

  # Remove legend if required
  if (legend == FALSE) {
    plot <- plot + guides(fill = "none")
  }

  return(plot)
}

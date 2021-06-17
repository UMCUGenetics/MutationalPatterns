#' Plot strand bias per base substitution type per group
#'
#' @param strand_bias data.frame, result from strand_bias function
#' @param colors Optional color vector with 6 values for plotting
#' @param sig_type The type of significance to be used. Possible values:
#'              * 'fdr' False discovery rate.
#'              A type of multiple testing correction.;
#'              * 'p' for regular p values.
#' @return Barplot
#'
#' @import ggplot2
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## following mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#'
#' tissue <- c(
#'   "colon", "colon", "colon",
#'   "intestine", "intestine", "intestine",
#'   "liver", "liver", "liver"
#' )
#'
#' ## Perform the strand bias test.
#' strand_counts <- strand_occurrences(mut_mat_s, by = tissue)
#' strand_bias <- strand_bias_test(strand_counts)
#'
#' ## Plot the strand bias.
#' plot_strand_bias(strand_bias)
#'
#' ## Use multiple (max 3) significance cutoffs.
#' ## This will vary the number of significance stars.
#' strand_bias_multistars <- strand_bias_test(strand_counts,
#'   p_cutoffs = c(0.05, 0.01, 0.005),
#'   fdr_cutoffs = c(0.1, 0.05, 0.01)
#' )
#' plot_strand_bias(strand_bias_multistars)
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurrences}},
#' \code{\link{strand_bias_test}}
#' \code{\link{plot_strand}}
#'
#' @export

plot_strand_bias <- function(strand_bias, colors = NA, sig_type = c("fdr", "p")) {
  

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  type <- significant <- strand_1 <- strand_2 <- log2_ratio <- NULL
  sig_plot <- log2_ratio_no1pseudo <- NULL

  # Match argument
  sig_type <- match.arg(sig_type)
  
  # if colors parameter not provided, set to default colors
  if (.is_na(colors)) {
    colors <- COLORS6
  }

  # get variable names
  var_names <- colnames(strand_bias)[3:4]
  colnames(strand_bias)[3:4] <- c("strand_1", "strand_2")


  # determine max y value for plotting
  # = log2 ratio with pseudo counts of 0.1
  strand_bias <- dplyr::mutate(strand_bias,
    log2_ratio = log2((strand_1 + 0.1) / (strand_2 + 0.1)),
    log2_ratio_no1pseudo = log2((strand_1) / (strand_2 + 0.1))
  )

  if (sig_type == "p") {
    strand_bias$sig_plot <- strand_bias$significant
  } else {
    strand_bias$sig_plot <- strand_bias$significant_fdr
  }

  # max yvalue for plotting plus
  max <- round(max(abs(strand_bias$log2_ratio)), digits = 1) + 0.1
  pos_stars <- abs(strand_bias$log2_ratio_no1pseudo)
  max_pos_star <- round(max(pos_stars[is.finite(pos_stars)]), digits = 1) + 0.1
  if (max < max_pos_star) {
    max <- max_pos_star
  }

  # add label for infinite values
  label2 <- log2(strand_bias$ratio)
  select <- which(is.finite(label2))
  label2[select] <- " "

  # plot strand bias with poisson test results
  plot <- ggplot(strand_bias, aes(x = type, y = log2_ratio, fill = type)) +
    geom_bar(colour = "black", stat = "identity", position = "identity") +
    scale_fill_manual(values = COLORS6) +
    scale_y_continuous(limits = c(-max, max)) +
    geom_text(
      aes(
        x = type,
        y = log2_ratio_no1pseudo,
        label = sig_plot,
        vjust = ifelse(sign(log2_ratio_no1pseudo) > 0, 0.5, 1)
      ),
      size = 8,
      position = ggplot2::position_dodge(width = 1)
    ) +
    facet_grid(. ~ group) +
    theme_bw() +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      legend.title = element_blank()
    ) +
    xlab("") +
    ylab(paste("log2(", var_names[1], "/", var_names[2], ")", sep = "")) +
    scale_x_discrete(breaks = NULL)

  return(plot)
}

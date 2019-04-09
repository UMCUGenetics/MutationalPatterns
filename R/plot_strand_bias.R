#' Plot strand bias per base substitution type per group
#'
#' @param strand_bias data.frame, result from strand_bias function
#' @param colors Optional color vector with 6 values for plotting
#' @return Barplot
#'
#' @import ggplot2
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## following mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' tissue <- c("colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver")
#'
#' ## Perform the strand bias test.
#' strand_counts = strand_occurrences(mut_mat_s, by=tissue)
#' strand_bias = strand_bias_test(strand_counts)
#'
#' ## Plot the strand bias.
#' plot_strand_bias(strand_bias)
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurrences}},
#' \code{\link{strand_bias_test}}
#' \code{\link{plot_strand}}
#'
#' @export

plot_strand_bias = function(strand_bias, colors)
{
  # get variable names
  var_names = colnames(strand_bias)[3:4]

  # if colors parameter not provided, set to default colors
  if (missing(colors))
    colors=COLORS6

  # determine max y value for plotting
  # = log2 ratio with pseudo counts of 0.1
  log2_ratio = log2(  (strand_bias[,3]+0.1) /
                        (strand_bias[,4]+0.1))

  # max yvalue for plotting plus
  max = round(max(abs(log2_ratio)), digits = 1) + 0.1
  pos_stars = abs(log2((strand_bias[, 3])/(strand_bias[, 4] + 0.1)))
  max_pos_star = round(max(pos_stars[is.finite(pos_stars)]), digits = 1) + 0.1
  if(max < max_pos_star){
    max = max_pos_star
}

  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL
  type = NULL
  significant = NULL

  # add label for infinite values
  label2 = log2(strand_bias$ratio)
  select = which(is.finite(label2))
  label2[select] = " "

  # plot strand bias with poisson test results
  plot = ggplot(strand_bias, aes( x = type,
                                  y = log2((strand_bias[,3]+0.1) /
                                             (strand_bias[,4]+0.1)),
                                  fill = type)) +
    scale_fill_manual(values = COLORS6) +
    geom_bar(colour = "black", stat ="identity", position = "identity") +
    scale_y_continuous(limits = c(-max, max)) +
    geom_text(
      aes(x = type,
          y = log2((strand_bias[,3]) / (strand_bias[,4]+0.1)),
          ymax = log2((strand_bias[,3]) / (strand_bias[,4]+0.1)),
          label = significant,
          vjust = ifelse(sign(log2((strand_bias[,3]) /
                                     (strand_bias[,4]+0.1))) > 0, 0.5, 1)),
      size = 8,
      position = ggplot2::position_dodge(width = 1)) +
    facet_grid(. ~ group) +
    theme_bw()  +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_blank()) +
    xlab("") +
    ylab(paste("log2(", var_names[1], "/", var_names[2], ")", sep = "") ) +
    scale_x_discrete(breaks=NULL)

    return(plot)
}

#' Plot strand per base substitution type
#'
#' For each base substitution type and transcriptional strand the total number
#' of mutations and the relative contribution within a group is returned.
#' @param strand_bias_df data.frame, result from strand_bias function
#' @param mode Either "absolute" for absolute number of mutations, or
#' "relative" for relative contribution, default = "relative"
#' @param colors Optional color vector for plotting with 6 values
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
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' tissue <- c(
#'   "colon", "colon", "colon",
#'   "intestine", "intestine", "intestine",
#'   "liver", "liver", "liver"
#' )
#'
#' strand_counts <- strand_occurrences(mut_mat_s, by = tissue)
#'
#' ## Plot the strand in relative mode.
#' strand_plot <- plot_strand(strand_counts)
#'
#' ## Or absolute mode.
#' strand_plot <- plot_strand(strand_counts, mode = "absolute")
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurrences}},
#' \code{\link{plot_strand_bias}}
#'
#' @export

plot_strand <- function(strand_bias_df, mode = c("relative", "absolute"), colors = NA) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  type <- relative_contribution <- no_mutations <- y_vals <- NULL

  # Match argument
  mode <- match.arg(mode)
  
  # if colors parameter not provided, set to default colors.
  if (.is_na(colors)) {
    colors <- COLORS6
  }
  
  # Set y values based on the plotting mode.
  if (mode == "relative") {
    strand_bias_df$y_vals <- strand_bias_df$relative_contribution
    y_lab <- "Relative contribution"
  } else if (mode == "absolute") {
    strand_bias_df$y_vals <- strand_bias_df$no_mutations
    y_lab <- "Total number of mutations"
  }

  # Create plot. Warnings about using alpha as a discrete variable are muffled.
  withCallingHandlers(
    {
      plot <- ggplot(strand_bias_df, aes(
        x = type,
        y = y_vals,
        fill = type,
        alpha = strand
      )) +
        geom_bar(
          stat = "identity",
          position = "dodge",
          colour = "black",
          cex = 0.5
        ) +
        scale_fill_manual(values = colors) +
        scale_alpha_discrete(range = c(1, 0.4)) +
        labs(y = y_lab, x = "") +
        facet_grid(. ~ group) +
        theme_bw() +
        scale_x_discrete(breaks = NULL)
    },
    warning = function(w) {
      if (grepl("Using alpha for a discrete variable is not advised.", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  return(plot)
}

#' Plot enrichment/depletion of mutations in genomic regions
#'
#' @param df Dataframe result from enrichment_depletion_test()
#' @param sig_type The type of significance to be used. Possible values:
#'              * 'fdr' False discovery rate.
#'              A type of multiple testing correction.;
#'              * 'p' for regular p values.
#' @return Plot with two parts. 1: Barplot with no. mutations expected and
#' observed per region. 2: Effect size of enrichment/depletion
#' (log2ratio) with results significance test.
#'
#' @import ggplot2
#'
#' @examples
#' ## See the 'genomic_distribution()' example for how we obtained the
#' ## following data:
#' distr <- readRDS(system.file("states/distr_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' tissue <- c(
#'   "colon", "colon", "colon",
#'   "intestine", "intestine", "intestine",
#'   "liver", "liver", "liver"
#' )
#'
#' ## Perform the enrichment/depletion test.
#' distr_test <- enrichment_depletion_test(distr, by = tissue)
#'
#' ## Plot the enrichment/depletion
#' plot_enrichment_depletion(distr_test)
#' 
#' #Perform and plot the enrichmet depletion test for all samples pooled
#' distr_test2 <- enrichment_depletion_test(distr)
#' plot_enrichment_depletion(distr_test2)
#'
#' ## Plot with p values instead of fdr
#' plot_enrichment_depletion(distr_test, sig_type = "p")
#'
#' ## Use multiple (max 3) significance cutoffs.
#' ## This will vary the number of significance stars.
#' distr_multistars <- enrichment_depletion_test(distr,
#'   by = tissue,
#'   p_cutoffs = c(0.05, 0.01, 0.005),
#'   fdr_cutoffs = c(0.1, 0.05, 0.01)
#' )
#' plot_enrichment_depletion(distr_multistars)
#' @seealso
#' \code{\link{enrichment_depletion_test}},
#' \code{\link{genomic_distribution}}
#'
#' @export

plot_enrichment_depletion <- function(df, sig_type = c("fdr", "p")) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  value <- variable <- observed <- expected <- significant <- NULL
  region <- ratio <- log2_ratio <- sig_plot <- NULL

  # Match argument
  sig_type <- match.arg(sig_type)

  # Create dataframe for absolute part of plot
  df_abs <- df %>%
    dplyr::select(by, region, observed, expected) %>%
    tidyr::pivot_longer(c(-by, -region), names_to = "variable", values_to = "value") %>%
    dplyr::mutate(variable = factor(variable, levels = unique(variable)))

  # Create dataframe for ratio part of plot
  df_ratio <- df %>%
    dplyr::mutate(
      ratio = (observed + 0.1) / (expected + 0.1),
      log2_ratio = log2(ratio)
    )

  if (sig_type == "p") {
    df_ratio$sig_plot <- df_ratio$significant
  } else {
    df_ratio$sig_plot <- df_ratio$significant_fdr
  }



  # Part 1: No. mutations expected and observed per region
  withCallingHandlers(
    {
      plot1 <- ggplot(df_abs, aes(
        x = by,
        y = value,
        fill = by,
        group = variable,
        alpha = variable
      )) +
        geom_bar(
          colour = "black",
          stat = "identity",
          position = position_dodge()
        ) +
        facet_grid(~region) +
        labs(x = "", y = "No. mutations") +
        scale_x_discrete(breaks = NULL) +
        scale_alpha_discrete(range = c(0.1, 1)) +
        theme_bw() +
        theme(
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_blank()
        )
    },
    warning = function(w) {
      if (grepl("Using alpha for a discrete variable is not advised.", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )

  # determine max y value for plotting
  # = log2 ratio with pseudo counts
  y_max <- round(max(abs(df_ratio$log2_ratio)), digits = 1) + 0.1

  # Part 2: effect size of enrichment/depletion with significance test
  plot2 <- ggplot(data = df_ratio, aes(
    x = by,
    y = log2_ratio,
    fill = by
  )) +
    geom_bar(
      colour = "black",
      stat = "identity",
      position = position_dodge()
    ) +
    scale_y_continuous(limits = c(-y_max, y_max)) +
    geom_text(
      aes(
        x = by,
        y = log2_ratio,
        label = sig_plot,
        vjust = ifelse(sign(log2_ratio) > 0, 0.5, 1)
      ),
      size = 8, position = position_dodge(width = 1)
    ) +
    facet_grid(~region) +
    labs(x = "", y = "log2(observed/expected)") +
    scale_x_discrete(breaks = NULL) +
    theme_bw() +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      legend.title = element_blank()
    )

  output <- cowplot::plot_grid(plot1, plot2, ncol = 1, nrow = 2, rel_heights = c(2, 1.2))
  return(output)
}

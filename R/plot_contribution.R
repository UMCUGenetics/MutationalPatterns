#' Plot signature contribution barplot
#'
#' Plot contribution of signatures. Can be used on both the results of a NMF
#' and on the results of signature refitting.
#'
#' @param contribution Signature contribution matrix
#' @param signatures Signature matrix.
#' Necessary when plotting NMF results in "absolute" mode.
#' It's not necessary in relative mode or when visualizing signature refitting results
#' @param index optional sample subset parameter
#' @param coord_flip Flip X and Y coordinates, default = FALSE
#' @param mode "relative" or "absolute"; to plot the relative contribution or
#' absolute number of mutations, default = "relative"
#' @param palette A color palette like c("#FF0000", "#00FF00", "9999CC") that
#' will be used as colors in the plot.  By default, ggplot2's colors are used
#' to generate a palette.
#'
#' @return Stacked barplot with contribution of each signature for each sample
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#' @examples
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Optionally set column and row names.
#' colnames(nmf_res$signatures) <- c("Signature A", "Signature B")
#' rownames(nmf_res$contribution) <- c("Signature A", "Signature B")
#'
#' ## Plot the relative contribution
#' plot_contribution(nmf_res$contribution)
#'
#' ## Plot the absolute contribution.
#' ## When plotting absolute NMF results, the signatures need to be included.
#' plot_contribution(nmf_res$contribution,
#'   nmf_res$signature,
#'   mode = "absolute"
#' )
#'
#'
#' ## Only plot a subset of samples
#' plot_contribution(nmf_res$contribution,
#'   nmf_res$signature,
#'   mode = "absolute",
#'   index = c(1, 2)
#' )
#' ## Flip the coordinates
#' plot_contribution(nmf_res$contribution,
#'   nmf_res$signature,
#'   mode = "absolute",
#'   coord_flip = TRUE
#' )
#'
#' ## You can also use the results of signature refitting.
#' ## Here we load some data as an example
#' fit_res <- readRDS(system.file("states/snv_refit.rds",
#'   package = "MutationalPatterns"
#' ))
#' plot_contribution(fit_res$contribution)
#'
#' ## Or again in absolute mode
#' plot_contribution(fit_res$contribution,
#'   mode = "absolute"
#' )
#' @seealso
#' \code{\link{extract_signatures}},
#' \code{\link{mut_matrix}}
#'
#' @export

plot_contribution <- function(contribution,
                              signatures = NA,
                              index = NA,
                              coord_flip = FALSE,
                              mode = c("relative", "absolute"),
                              palette = NA) {
  # Match argument
  mode <- match.arg(mode)

  # optional subsetting if index parameter is provided
  if (!.is_na(index)) {
    contribution <- contribution[, index, drop = FALSE]
  }

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  Sample <- Contribution <- Signature <- NULL

  # When working on NMF results, the contribution needs to be multiplied by the signature colSums.
  if (mode == "absolute" & !.is_na(signatures)) {
    # calculate signature contribution in absolute number of signatures
    total_signatures <- colSums(signatures)
    abs_contribution <- contribution * total_signatures
  }

  # Make data long. Also create factors for ordering.
  tb <- contribution %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Signature") %>%
    tidyr::pivot_longer(-Signature, names_to = "Sample", values_to = "Contribution") %>%
    dplyr::mutate(
      Sample = factor(Sample, levels = unique(Sample)),
      Signature = factor(Signature, levels = unique(Signature))
    )

  # Different plotting between absolute and relative
  if (mode == "absolute") {
    bar_geom <- geom_bar(stat = "identity", colour = "black")
    y_lab <- "Absolute contribution \n (no. mutations)"
  } else if (mode == "relative") {
    bar_geom <- geom_bar(position = "fill", stat = "identity", colour = "black")
    y_lab <- "Relative contribution"
  }
  
  # Determine what signatures are present for the legend.
  present_sigs <- tb %>% 
    dplyr::filter(Contribution != 0) %>% 
    dplyr::pull(Signature) %>% 
    unique()
  
  #Create plot
  plot <- ggplot(tb, aes(x = Sample, y = Contribution, fill = Signature)) +
    bar_geom +
    labs(x = "", y = y_lab) +
    scale_fill_discrete(breaks = present_sigs) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
    )

  # Allow custom color palettes.
  if (!.is_na(palette)) {
    plot <- plot + scale_fill_manual(name = "Signature", values = palette)
  }

  # Handle coord_flip.
  if (coord_flip) {
    plot <- plot + coord_flip() + xlim(rev(levels(factor(tb$Sample))))
  } else {
    plot <- plot + xlim(levels(factor(tb$Sample)))
  }

  return(plot)
}

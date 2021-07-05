#' Plot 96 trinucleotide profile
#'
#' Plot relative contribution of 96 trinucleotides
#' @param mut_matrix 96 trinucleotide profile matrix
#' @param ymax Y axis maximum value, default = 0.2
#' @param colors Optional 6 value color vector.
#' @param condensed More condensed plotting format. Default = F.
#' @return 96 trinucleotide profile plot
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the
#' ## mutation matrix information:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Plot the 96-profile of three samples
#' plot_96_profile(mut_mat[, c(1, 4, 7)])
#'
#' ## Plot a condensed profile
#' plot_96_profile(mut_mat[, c(1, 4, 7)], condensed = TRUE)
#'
#' ## It's also possible to plot signatures, for example signatures
#' ## generated with NMF
#' ## See 'extract_signatures()' on how we obtained these signatures.
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Optionally, provide signature names
#' colnames(nmf_res$signatures) <- c("Signature A", "Signature B")
#'
#' ## Generate the plot
#' plot_96_profile(nmf_res$signatures)
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{plot_profile_heatmap}},
#' \code{\link{plot_river}}
#'
#' @export

plot_96_profile <- function(mut_matrix, colors = NA, ymax = 0.2, condensed = FALSE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  freq <- full_context <- substitution <- context <- NULL

  # Check color vector length
  # Colors for plotting
  if (.is_na(colors)) {
    colors <- COLORS6
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }

  # Make contribution relative
  norm_mut_matrix <- apply(mut_matrix, 2, function(x) x / sum(x))

  # Get substitution and context from rownames and make long.
  tb <- norm_mut_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("full_context") %>%
    dplyr::mutate(
      substitution = stringr::str_replace(full_context, "\\w\\[(.*)\\]\\w", "\\1"),
      context = stringr::str_replace(full_context, "\\[.*\\]", "\\.")
    ) %>%
    dplyr::select(-full_context) %>%
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", values_to = "freq") %>% 
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))

  # Change plotting parameters based on whether plot should be condensed.
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  } else {
    width <- 0.6
    spacing <- 0.5
  }

  # Create figure
  plot <- ggplot(data = tb, aes(
    x = context,
    y = freq,
    fill = substitution,
    width = width
  )) +
    geom_bar(stat = "identity", colour = "black", size = .2) +
    scale_fill_manual(values = colors) +
    facet_grid(sample ~ substitution) +
    ylab("Relative contribution") +
    coord_cartesian(ylim = c(0, ymax)) +
    scale_y_continuous(breaks = seq(0, ymax, 0.1)) +
    guides(fill = "none") +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 12, vjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = unit(spacing, "lines")
    )

  return(plot)
}

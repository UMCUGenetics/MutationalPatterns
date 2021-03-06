#' Plot the DBS contexts
#'
#' @details
#' Plots the number of DBS COSMIC context per sample.
#' It takes a tibble with counts as its input. This tibble can be generated by count_dbs_contexts
#' Each sample is plotted in a separate facet.
#' The same y axis can be used for all samples or a separate y axis can be used.
#'
#' @param counts A tibble containing the number of DBS per COSMIC context.
#' @param same_y A boolean describing whether the same y axis should be used for all samples.
#' @param condensed More condensed plotting format. Default = F.
#'
#' @return A ggplot figure.
#'
#' @examples
#' ## Get The DBS counts
#' ## See 'count_dbs_contexts()' for more info on how to do this.
#' dbs_counts <- readRDS(system.file("states/blood_dbs_counts.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Plot contexts
#' plot_dbs_contexts(dbs_counts)
#'
#' ## Use the same y axis for all samples.
#' plot_dbs_contexts(dbs_counts, same_y = TRUE)
#'
#' ## Create a more condensed plot
#' plot_dbs_contexts(dbs_counts, condensed = TRUE)
#' @import ggplot2
#' @importFrom magrittr %>%
#' @family DBS
#'
#' @seealso \code{\link{count_dbs_contexts}}, \code{\link{plot_main_dbs_contexts}}
#'
#' @export
plot_dbs_contexts <- function(counts, same_y = FALSE, condensed = FALSE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  count <- REF <- ALT <- muttype_total <- sample <- NULL

  # Transform to data frame
  counts <- counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("muttype_total") %>%
    tidyr::separate(muttype_total, c("REF", "ALT"), sep = "_") %>%
    dplyr::mutate(REF = factor(REF, levels = BiocGenerics::unique(REF)))

  # Set levels of ALT
  bases <- c("A", "C", "G", "T")
  bases1 <- bases
  bases_combi <- tidyr::crossing(bases, bases1)
  counts$ALT <- factor(counts$ALT, levels = stringr::str_c(bases_combi$bases, bases_combi$bases1))


  # Transform data to long format.
  counts <- counts %>% 
    tidyr::gather(key = "sample", value = "count", -REF, -ALT) %>% 
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))

  # Count nr of mutations
  nr_muts <- counts %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(nr_muts = round(sum(count)))

  if (same_y) {
    facet_scale <- "free_x"
  } else {
    facet_scale <- "free"
  }

  # Create facet labs
  facet_labs_y <- stringr::str_c(nr_muts$sample, " (n = ", nr_muts$nr_muts, ")")
  names(facet_labs_y) <- nr_muts$sample
  facet_labs_x <- stringr::str_c(levels(counts$REF), ">NN")
  names(facet_labs_x) <- levels(counts$REF)

  # Change plotting parameters based on whether plot should be condensed.
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  } else {
    width <- 0.6
    spacing <- 0.5
  }

  # Set colours
  colors <- c(
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"
  )

  # Create plot
  fig <- ggplot(counts, aes(x = ALT, y = count, fill = REF, width = width)) +
    geom_bar(stat = "identity") +
    facet_grid(sample ~ REF,
      scales = facet_scale,
      space = "free_x",
      labeller = labeller(REF = facet_labs_x, sample = facet_labs_y)
    ) +
    scale_fill_manual(guide = FALSE, values = colors) +
    labs(fill = "Mutation type", title = "", y = "Nr of DBSs", x = "") +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing.x = unit(spacing, "lines"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  return(fig)
}

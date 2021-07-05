#' Compare two DBS mutation profiles
#'
#' Plots two DBS mutation profiles and their difference, reports the residual
#' sum of squares (RSS).
#'
#' @param profile1 First mutation profile
#' @param profile2 Second mutation profile
#' @param profile_names Character vector with names of the mutations profiles
#' used for plotting, default = c("profile 1", "profile 2")
#' @param profile_ymax Maximum value of y-axis (relative contribution) for
#' profile plotting. This can only be used to increase the y axis.
#' If bars fall outside this limit, the maximum value is
#' automatically increased. default = 0.2.
#' @param diff_ylim Y-axis limits for profile difference plot,
#' default = c(-0.1, 0.1)
#'
#' @return A ggplot2 object
#' @export
#' @family DBS
#' @seealso \code{\link{plot_compare_profiles}},
#' \code{\link{plot_compare_indels}},
#' \code{\link{plot_compare_mbs}}
#'
#' @examples
#'
#' ## Get the DBS counts
#' ## See 'count_dbs_contexts()' for more info on how to do this.
#' dbs_counts <- readRDS(system.file("states/blood_dbs_counts.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Get DBS refit info.
#' ## See 'fit_to_signatures()' for more info on how to do this.
#' fit_res <- readRDS(system.file("states/dbs_refit.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Compare the reconstructed profile of sample 1 with the original profile
#' ## The same thing could be done with a reconstructed profile from NMF.
#' plot_compare_dbs(dbs_counts[, 1], fit_res$reconstructed[, 1])
#'
#' ## You could also compare regular mutation profiles with eachother.
#' plot_compare_dbs(
#'   dbs_counts[, 1],
#'   dbs_counts[, 2]
#' )
#'
#' ## Or change the names of the profiles
#' plot_compare_dbs(dbs_counts[, 1],
#'   dbs_counts[, 2],
#'   profile_names = c("Original", "Reconstructed")
#' )
#'
#' ## You can also change the y limits.
#' ## This can be done separately for the profiles and the different facets.
#' plot_compare_dbs(dbs_counts[, 1],
#'   dbs_counts[, 2],
#'   profile_ymax = 0.3,
#'   diff_ylim = c(-0.03, 0.03)
#' )
plot_compare_dbs <- function(profile1, profile2,
                             profile_names = c("profile 1", "profile 2"),
                             profile_ymax = 0.2,
                             diff_ylim = c(-0.1, 0.1)) {
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  count <- REF <- ALT <- value <- muttype_total <- NULL


  # Create a comparison of the profiles.
  comp <- .create_profile_comparison(profile1, profile2, profile_names)

  # Transform to data frame
  counts <- comp$matrix %>%
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
  counts <- tidyr::gather(counts, key = "sample", value = "count", -REF, -ALT) %>%
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))

  # Add dummy non_visible data points to force y axis limits per facet
  df_blank <- .create_dummy_limits(counts[, c("REF", "ALT")], profile_names, profile_ymax, diff_ylim)

  # Set facet labels
  facet_labs_x <- stringr::str_c(levels(counts$REF), ">NN")
  names(facet_labs_x) <- levels(counts$REF)

  # Create plot
  fig <- ggplot(counts, aes(x = ALT, y = count, fill = REF)) +
    geom_bar(stat = "identity") +
    geom_blank(data = df_blank, aes(x = ALT, y = value)) +
    facet_grid(sample ~ REF,
      scales = "free",
      space = "free_x",
      labeller = labeller(REF = facet_labs_x)
    ) +
    scale_fill_manual(guide = "none", values = DBS_COLORS) +
    labs(fill = "Mutation type", title = comp$title, y = "Relative contribution", x = "") +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "cadetblue"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  return(fig)
}

#' Compare two mbs mutation profiles
#'
#' Plots two mbs mutation profiles and their difference, reports the residual
#' sum of squares (RSS).
#'
#' @param profile1 First mutation profile
#' @param profile2 Second mutation profile
#' @param profile_names Character vector with names of the mutations profiles
#' used for plotting, default = c("profile 1", "profile 2")
#' @param profile_ymax Maximum value of y-axis (relative contribution) for
#' profile plotting. This can only be used to increase the y axis.
#' If bars fall outside this limit, the maximum value is
#' automatically increased. default = 1.
#' @param diff_ylim Y-axis limits for profile difference plot,
#' default = c(-0.5, 0.5)
#'
#' @return A ggplot2 object
#' @export
#' @family MBS
#' @seealso \code{\link{plot_compare_profiles}},
#' \code{\link{plot_compare_dbs}},
#' \code{\link{plot_compare_indels}}
#' @examples
#'
#' ## Get the mbs counts
#' ## See 'count_mbs_contexts()' for more info on how to do this.
#' mbs_counts <- readRDS(system.file("states/blood_mbs_counts.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#'
#' ## You could compare regular mutation profiles with eachother.
#' plot_compare_mbs(
#'   mbs_counts[, 1],
#'   mbs_counts[, 2]
#' )
#'
#' ## Or change the names of the profiles
#' plot_compare_mbs(mbs_counts[, 1],
#'   mbs_counts[, 2],
#'   profile_names = c("Original", "Reconstructed")
#' )
#'
#' ## You can also change the y limits.
#' ## This can be done separately for the profiles and the different facets.
#' plot_compare_mbs(mbs_counts[, 1],
#'   mbs_counts[, 2],
#'   profile_ymax = 0.9,
#'   diff_ylim = c(-0.8, 0.8)
#' )
#'
#' ## You could also compare a reconstructed profile.
#' ## However, the example data does not contain enough MBS variants to use NMF.
#' ## Existing signatures have also not yet been defined.
plot_compare_mbs <- function(profile1, profile2,
                             profile_names = c("profile 1", "profile 2"),
                             profile_ymax = 1,
                             diff_ylim = c(-0.5, 0.5)) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  count <- size <- value <- muttype_total <- NULL

  # Create a comparison of the profiles.
  comp <- .create_profile_comparison(profile1, profile2, profile_names)

  # Transform to data frame
  counts <- comp$matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("size") %>%
    tidyr::pivot_longer(-size, names_to = "sample", values_to = "count") %>%
    dplyr::mutate(
      size = factor(size, levels = unique(size)),
      sample = factor(sample, levels = unique(sample))
    )

  # Add dummy non_visible data points to force y axis limits per facet
  df_blank <- .create_dummy_limits(counts[, c("size")], profile_names, profile_ymax, diff_ylim)

  # Create plot
  fig <- ggplot(counts, aes(x = size, y = count, fill = size)) +
    geom_bar(stat = "identity") +
    geom_blank(data = df_blank, aes(x = size, y = value)) +
    facet_grid(sample ~ .,
      scales = "free_y"
    ) +
    labs(x = "MBS size", y = "Relative contribution", title = comp$title) +
    scale_fill_manual(values = MBS_COLORS, guide = "none") +
    theme_classic() +
    theme(
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA),
      axis.text = element_text(
        size = rel(0.8),
        colour = "black"
      ),
      axis.ticks = element_line(colour = "black", size = 1),
      axis.line = element_line(size = 1),
      strip.background = element_blank()
    )
  return(fig)
}

#' Compare two 96 mutation profiles
#'
#' Plots two 96 mutation profiles and their difference, reports the residual
#' sum of squares (RSS).
#'
#' @param profile1 First 96 mutation profile
#' @param profile2 Second 96 mutation profile
#' @param profile_names Character vector with names of the mutations profiles
#' used for plotting, default = c("profile 1", "profile 2")
#' @param profile_ymax Maximum value of y-axis (relative contribution) for
#' profile plotting. This can only be used to increase the y axis.
#' If bars fall outside this limit, the maximum value is
#' automatically increased. default = 0.2.
#' @param diff_ylim Y-axis limits for profile difference plot,
#' default = c(-0.02, 0.02)
#' @param colors 6 value color vector
#' @param condensed More condensed plotting format. Default = F.
#' @return 96 spectrum plot of profile 1, profile 2 and their difference
#'
#' @import ggplot2
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the following
#' ## mutation matrix.
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Compare the reconstructed 96-profile of sample 1 with the original profile
#' ## The same thing could be done with a reconstructed profile from signature refitting.
#' plot_compare_profiles(mut_mat[, 1],
#'   nmf_res$reconstructed[, 1],
#'   profile_names = c("Original", "Reconstructed")
#' )
#'
#' ## You could also compare regular mutation profiles with eachother.
#' plot_compare_profiles(
#'   mut_mat[, 1],
#'   mut_mat[, 2]
#' )
#'
#'
#' ## You can also change the y limits.
#' ## This can be done separately for the profiles and the different facets.
#' plot_compare_profiles(mut_mat[, 1],
#'   mut_mat[, 2],
#'   profile_ymax = 0.3,
#'   diff_ylim = c(-0.03, 0.03)
#' )
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{extract_signatures}},
#' \code{\link{plot_compare_indels}},
#' \code{\link{plot_compare_dbs}},
#' \code{\link{plot_compare_mbs}}
#'
#' @export
#'
plot_compare_profiles <- function(profile1,
                                  profile2,
                                  profile_names = c("profile 1", "profile 2"),
                                  profile_ymax = 0.2,
                                  diff_ylim = c(-0.02, 0.02),
                                  colors = NA,
                                  condensed = FALSE) {
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  value <- substitution <- Sample <- Contribution <- Signature <- variable <- NULL
  full_context <- context <- NULL

  # if colors parameter not provided, set to default colors
  if (.is_na(colors)) {
    colors <- COLORS6
  }

  # Create a comparison of the profiles.
  comp <- .create_profile_comparison(profile1, profile2, profile_names)


  # Get substitution and context. Then make data long for plotting.
  df <- comp$matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("full_context") %>%
    dplyr::mutate(
      substitution = stringr::str_replace(full_context, "\\w\\[(.*)\\]\\w", "\\1"),
      context = stringr::str_replace(full_context, "\\[.*\\]", "\\.")
    ) %>%
    dplyr::select(-full_context) %>%
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", values_to = "value") %>%
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))


  # Add dummy non_visible data points to force y axis limits per facet
  df_blank <- .create_dummy_limits(df[, c("substitution", "context")], profile_names, profile_ymax, diff_ylim)

  # Plotting parameters
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  } else {
    width <- 0.6
    spacing <- 0.5
  }

  # Create plot
  plot <- ggplot(data = df, aes(
    x = context,
    y = value,
    fill = substitution,
    width = width
  )) +
    geom_bar(
      stat = "identity",
      position = "identity",
      colour = "black", size = .2
    ) +
    geom_blank(data = df_blank, aes(x = context, y = value)) +
    scale_fill_manual(values = colors) +
    facet_grid(sample ~ substitution, scales = "free_y") +
    labs(
      y = "Relative contribution",
      title = comp$title
    ) +
    guides(fill = "none") +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 12, vjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = unit(spacing, "lines")
    )

  return(plot)
}


#' Create a relative comparison between two profiles.
#'
#' Create a matrix with the relative profiles and the difference.
#'
#' @param profile1 First mutation profile
#' @param profile2 Second mutation profile
#' @param profile_names Character vector with names of the mutations profiles
#' used for plotting
#'
#' @return matrix with the relative profiles and the difference
#' @noRd
#'
.create_profile_comparison <- function(profile1, profile2, profile_names) {
  s1_relative <- profile1 / sum(profile1)
  s2_relative <- profile2 / sum(profile2)
  diff <- s1_relative - s2_relative

  # residual sum of squares
  RSS <- sum(diff^2)
  RSS <- format(RSS, scientific = TRUE, digits = 3)

  # calculate cosine similarity between the two profiles
  cosine_sim <- cos_sim(profile1, profile2)
  # round
  cosine_sim <- round(cosine_sim, 3)

  # Create title
  title <- paste0("RSS = ", RSS, "; Cosine similarity = ", cosine_sim)

  # Combine samples and diff.
  x <- BiocGenerics::cbind(s1_relative, s2_relative, diff)
  colnames(x) <- c(profile_names, "Difference")

  res <- list("title" = title, "matrix" = x)
  return(res)
}


#' Create a dummy data frame with y-axis limits
#'
#' This functions creates dummy dataframe with y-axis limits for the relative profiles,
#' as well as the diff profile. The result can be used by ggplot to change the y axis
#' separately for the profiles and diffs.
#'
#' @param df Dataframe with the mutation types
#' @param profile_names Character vector with names of the mutations profiles
#' used for plotting
#' @param profile_ymax Maximum value of y-axis (relative contribution) for
#' profile plotting
#' @param diff_ylim Y-axis limits for profile difference plot
#'
#' @return Dataframe with y-axis limits
#' @noRd
#' @importFrom magrittr %>%
#'
.create_dummy_limits <- function(df, profile_names, profile_ymax, diff_ylim) {
  df_dummy <- df[c(1, 1, 1, 1), ] %>%
    dplyr::mutate(
      sample = c(profile_names, "Difference", "Difference"),
      sample = factor(sample, levels = unique(sample)),
      value = c(profile_ymax, profile_ymax, diff_ylim[1], diff_ylim[2])
    )
  return(df_dummy)
}

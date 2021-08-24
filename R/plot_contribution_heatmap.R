#' Plot signature contribution heatmap
#'
#' Plot relative contribution of signatures in a heatmap
#'
#' @param contribution Signature contribution matrix
#' @param sig_order Character vector with the desired order of the signature names for plotting. Optional.
#' @param sample_order Character vector with the desired order of the sample names for plotting. Optional.
#' @param cluster_samples Hierarchically cluster samples based on euclidean distance. Default = T.
#' @param cluster_sigs Hierarchically cluster sigs based on euclidean distance. Default = T.
#' @param method The agglomeration method to be used for hierarchical clustering. This should be one of
#' "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC)
#' or "centroid" (= UPGMC). Default = "complete".
#' @param plot_values Plot relative contribution values in heatmap. Default = F.
#'
#' @return Heatmap with relative contribution of each signature for each sample
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#' @examples
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Set signature names as row names in the contribution matrix
#' rownames(nmf_res$contribution) <- c("Signature A", "Signature B")
#'
#' ## Plot with clustering.
#' plot_contribution_heatmap(nmf_res$contribution, cluster_samples = TRUE, cluster_sigs = TRUE)
#'
#' ## Define signature and sample order for plotting. If you have a mutation or signature
#' ## matrix, then this can be done like in the example of 'plot_cosine_heatmap()'
#' sig_order <- c("Signature B", "Signature A")
#' sample_order <- c(
#'   "colon1", "colon2", "colon3", "intestine1", "intestine2",
#'   "intestine3", "liver3", "liver2", "liver1"
#' )
#' plot_contribution_heatmap(nmf_res$contribution,
#'   cluster_samples = FALSE,
#'   sig_order = sig_order, sample_order = sample_order
#' )
#'
#' ## It's also possible to create a contribution heatmap with text values
#' output_text <- plot_contribution_heatmap(nmf_res$contribution, plot_values = TRUE)
#'
#' ## This function can also be used on the result of a signature refitting analysis.
#' ## Here we load a existing result as an example.
#' snv_refit <- readRDS(system.file("states/strict_snv_refit.rds",
#'   package = "MutationalPatterns"
#' ))
#' plot_contribution_heatmap(snv_refit$contribution, cluster_samples = TRUE, cluster_sigs = TRUE)
#' @seealso
#' \code{\link{extract_signatures}},
#' \code{\link{mut_matrix}},
#' \code{\link{plot_contribution}},
#' \code{\link{plot_cosine_heatmap}}
#'
#' @export

# plotting function for relative contribution of signatures in heatmap
plot_contribution_heatmap <- function(contribution, sig_order = NA, sample_order = NA, cluster_samples = TRUE,
                                      cluster_sigs = FALSE, method = "complete", plot_values = FALSE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  Signature <- Sample <- Contribution <- x <- y <- xend <- yend <- NULL

  # check contribution argument
  if (!inherits(contribution, "matrix")) {
    stop("contribution must be a matrix")
  }
  # check if there are signatures names in the contribution matrix
  if (is.null(row.names(contribution))) {
    stop("contribution must have row.names (signature names)")
  }

  # transpose
  contribution <- t(contribution)
  # relative contribution
  contribution_norm <- contribution / rowSums(contribution)

  # If cluster_samples is TRUE perform clustering. Else use supplied sample_order or
  # the current column order.
  if (!.is_na(sample_order) & cluster_samples == TRUE) {
    stop("sample_order can only be provided when cluster_samples is FALSE", call. = FALSE)
  } else if (!.is_na(sample_order)) {
    # check sample_order argument
    if (!inherits(sample_order, "character")) {
      stop("sample_order must be a character vector", call. = FALSE)
    }
    if (length(sample_order) != nrow(contribution_norm)) {
      stop("sample_order must have the same length as the number
          of samples in the explained matrix", call. = FALSE)
    }
  } else if (cluster_samples == TRUE) {
    # cluster samples based on eucledian distance between relative contribution_norm
    hc.sample <- hclust(dist(contribution_norm), method = method)
    # order samples according to clustering
    sample_order <- rownames(contribution_norm)[hc.sample$order]

    dhc <- as.dendrogram(hc.sample)
    # rectangular lines
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram_rows <- ggplot(ggdendro::segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +
      scale_y_reverse(expand = c(0.2, 0)) +
      ggdendro::theme_dendro()
  }
  else {
    sample_order <- rownames(contribution_norm)
  }


  # If cluster_sigs is TRUE perform clustering. Else use supplied sig_order or
  # the current column order.
  if (!.is_na(sig_order) & cluster_sigs == TRUE) {
    stop("sig_order can only be provided when cluster_sigs is FALSE", call. = FALSE)
  } else if (!.is_na(sig_order)) {
    # check sig_order argument
    if (!inherits(sig_order, "character")) {
      stop("sig_order must be a character vector", call. = FALSE)
    }
    if (length(sig_order) != ncol(contribution_norm)) {
      stop("sig_order must have the same length as the number
           of signatures in the explained matrix", call. = FALSE)
    }
  } else if (cluster_sigs == TRUE) {
    # Cluster cols
    hc.sample2 <- contribution_norm %>%
      t() %>%
      dist() %>%
      hclust(method = method)
    sig_order <- colnames(contribution_norm)[hc.sample2$order]

    dhc <- as.dendrogram(hc.sample2)
    # rectangular lines
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram_cols <- ggplot(ggdendro::segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      ggdendro::theme_dendro() +
      scale_y_continuous(expand = c(0.2, 0))
  } else {
    sig_order <- colnames(contribution_norm)
  }

  # Make matrix long and set factor levels, to get the correct order for plotting.
  contribution_norm.m <- contribution_norm %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(-Sample, names_to = "Signature", values_to = "Contribution") %>%
    dplyr::mutate(
      Signature = factor(Signature, levels = sig_order),
      Sample = factor(Sample, levels = sample_order)
    )

  # plot heatmap
  heatmap <- ggplot(contribution_norm.m, aes(x = Signature, y = Sample, fill = Contribution, order = Sample)) +
    geom_raster() +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Relative \ncontribution", limits = c(0, 1)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = NULL, y = NULL)
  # if plot_values is TRUE, add values to heatmap
  if (plot_values) {
    heatmap <- heatmap + geom_text(aes(label = round(Contribution, 2)), size = 3)
  }

  # if cluster_samples is TRUE, make dendrogram
  if (cluster_samples == TRUE & cluster_sigs == TRUE) {
    empty_fig <- ggplot() +
      theme_void()
    plot_final <- cowplot::plot_grid(empty_fig, dendrogram_cols, dendrogram_rows, heatmap,
      align = "hv", axis = "tblr", rel_widths = c(0.3, 1), rel_heights = c(0.3, 1)
    )
  }
  else if (cluster_samples == TRUE & cluster_sigs == FALSE) {
    # combine plots
    plot_final <- cowplot::plot_grid(dendrogram_rows, heatmap, align = "h", rel_widths = c(0.3, 1))
  } else if (cluster_samples == FALSE & cluster_sigs == TRUE) {
    plot_final <- cowplot::plot_grid(dendrogram_cols, heatmap, align = "v", rel_heights = c(0.3, 1)) +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(contribution_norm.m$Sample))))
  } else {
    plot_final <- heatmap +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(contribution_norm.m$Sample))))
  }

  return(plot_final)
}

#' Plot cosine similarity heatmap
#' 
#' Plot pairwise cosine similarities in a heatmap.
#' 
#' 
#' @param cos_sim_matrix Matrix with pairwise cosine similarities.
#'                       Result from \code{\link{cos_sim_matrix}}
#' @param col_order Character vector with the desired order of the columns names for plotting. Optional.
#' @param cluster_rows Hierarchically cluster rows based on eucledian distance. Default = TRUE.
#' @param method The agglomeration method to be used for hierarchical clustering. This should be one of 
#' "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) 
#' or "centroid" (= UPGMC). Default = "complete".
#' @param plot_values Plot cosine similarity values in heatmap. Default = FALSE.
#'
#' @return Heatmap with cosine similarities
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggdendro dendro_data segment theme_dendro
#' @importFrom cowplot plot_grid
#'
#' @examples
#' 
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## You can download the signatures from the COSMIC:
#' # http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#' 
#' ## Match the order to MutationalPatterns standard of mutation matrix
#' order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
#' ## Reorder cancer signatures dataframe
#' cancer_signatures = cancer_signatures[order,]
#' ## Use trinucletiode changes names as row.names
#' ## row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
#' ## Keep only 96 contributions of the signatures in matrix
#' cancer_signatures = as.matrix(cancer_signatures[,4:33])
#' ## Rename signatures to number only
#' colnames(cancer_signatures) = as.character(1:30)
#'
#' ## Calculate the cosine similarity between each signature and each 96 mutational profile
#' cos_matrix = cos_sim_matrix(mut_mat, cancer_signatures)
#' 
#' ## Cluster signatures based on cosine similarity 
#' sig_hclust = cluster_signatures(cancer_signatures)
#' col_order = colnames(cancer_signatures)[sig_hclust$order]
#' 
#' ## Plot the cosine similarity between each signature and each sample with hierarchical 
#' ## sample clustering and signatures order based on similarity
#' 
#' plot_cosine_heatmap(cos_matrix, col_order, cluster_rows = TRUE, method = "complete")
#' 
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{cos_sim_matrix}},
#' \code{\link{cluster_signatures}} 
#' 
#' @export

plot_cosine_heatmap = function(cos_sim_matrix, col_order, cluster_rows = TRUE, method = "complete", plot_values = FALSE)
{
  # check explained argument
  if(class(cos_sim_matrix) != "matrix")
  {stop("cos_sim_matrix must be a matrix")}
  # matrix should have row and colnames
  if(length(colnames(cos_sim_matrix)) == 0)
  {stop("cos_sim_matrix is missing colnames")}
  if(length(rownames(cos_sim_matrix)) == 0)
  {stop("cos_sim_matrix is missing rownames")}
  # if no signature order is provided, use the order as in the input matrix
  if(missing(col_order))
  {
    col_order = colnames(cos_sim_matrix)
  }
  # check col_order argument
  if(class(col_order) != "character")
  {stop("col_order must be a character vector")}
  if(length(col_order) != ncol(cos_sim_matrix))
  {stop("col_order must have the same length as the number of signatures in the explained matrix")}
  
  # if cluster samples is TRUE, perform clustering
  if(cluster_rows == TRUE)
  {
  # cluster samples based on eucledian distance between relative contribution
  hc.sample = hclust(dist(cos_sim_matrix), method = method)
  # order samples according to clustering
  sample_order = rownames(cos_sim_matrix)[hc.sample$order]
  }
  else
  {
    sample_order = rownames(cos_sim_matrix)
  }
  
  Cosine.sim = NULL
  Signature = NULL
  Sample = NULL
  x = NULL
  y = NULL
  xend = NULL
  yend = NULL

  # melt
  cos_sim_matrix.m = melt(cos_sim_matrix)
  # assign variable names
  colnames(cos_sim_matrix.m) = c("Sample", "Signature", "Cosine.sim")
  
  # change factor levels to the correct order for plotting
  cos_sim_matrix.m$Signature = factor(cos_sim_matrix.m$Signature, levels = col_order)
  cos_sim_matrix.m$Sample = factor(cos_sim_matrix.m$Sample, levels = sample_order)
  # plot heatmap
  heatmap = ggplot(cos_sim_matrix.m, aes(x=Signature, y=Sample, fill=Cosine.sim, order=Sample)) + 
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine \nsimilarity", limits = c(0,1)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x=NULL, y=NULL)
  # if plot_values is TRUE, add values to heatmap
  if (plot_values)
  {
    heatmap = heatmap + geom_text(aes(label = round(Cosine.sim, 2)), size = 3)
  }
  
  # if cluster samples is TRUE, make dendrogram
  if(cluster_rows == TRUE)
  {
    # get dendrogram
    dhc = as.dendrogram(hc.sample)
    # rectangular lines
    ddata = dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram = ggplot(segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + 
      scale_y_reverse(expand = c(0.2, 0)) + 
      theme_dendro()
    # combine plots
    plot_final = plot_grid(dendrogram, heatmap, align='h', rel_widths=c(0.3,1))
  }
  else
  {
    plot_final = heatmap +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  }
  
  return(plot_final)
}

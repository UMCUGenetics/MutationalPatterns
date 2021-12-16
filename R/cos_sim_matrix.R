#' Compute all pairwise cosine similarities between mutational profiles/signatures
#'
#' Computes all pairwise cosine similarities between the mutational profiles provided in the two mutation count matrices.
#' The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much two vectors are alike.
#'
#' @param mut_matrix1 mutation count matrix (dimensions: a mutation features X n samples)
#' @param mut_matrix2 96 mutation count matrix (dimensions: a mutation features X m samples)
#' @return Matrix with pairwise cosine similarities (dimensions: n mutational profiles X m mutational profiles)
#'
#' @examples
#' ## Get signatures
#' signatures <- get_known_signatures()
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#'
#' ## Calculate the cosine similarity between each COSMIC signature and each 96 mutational profile
#' cos_sim_matrix(mut_mat, signatures)
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{fit_to_signatures}},
#' \code{\link{plot_cosine_heatmap}}
#'
#' @export

cos_sim_matrix <- function(mut_matrix1, mut_matrix2) {
  
  # Check that both inputs are numeric.
  if (!all(apply(mut_matrix1, 2, is.numeric))){
    stop("The first input contains non-numeric columns, while all columns should be numeric.")
  }
  if (!all(apply(mut_matrix2, 2, is.numeric))){
    stop("The second input contains non-numeric columns, while all columns should be numeric.")
  }
  
  # Determine number of samples
  n_samples1 <- ncol(mut_matrix1)
  n_samples2 <- ncol(mut_matrix2)
  res_matrix <- matrix(nrow = n_samples1, ncol = n_samples2)

  # Loop over the columns of both input matrices,
  # to determine the cosine similarities.
  for (s in seq_len(n_samples1))
  {
    signal1 <- mut_matrix1[, s, drop = TRUE]
    cos_sim_vector <- c()
    for (i in seq_len(n_samples2))
    {
      signal2 <- mut_matrix2[, i, drop = TRUE]
      cos_sim_vector[i] <- cos_sim(signal1, signal2)
    }
    res_matrix[s, ] <- cos_sim_vector
  }
  rownames(res_matrix) <- colnames(mut_matrix1)
  colnames(res_matrix) <- colnames(mut_matrix2)

  return(res_matrix)
}

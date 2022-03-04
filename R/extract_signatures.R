#' Extract mutational signatures from 96 mutation matrix using NMF
#'
#' Decomposes trinucleotide count matrix into signatures and contribution of
#' those signatures to the spectra of the samples/vcf files.
#'
#' @param mut_matrix 96 mutation count matrix
#' @param rank Number of signatures to extract
#' @param nrun Number of iterations, default = 200.
#' A lower number will be faster, but result in less accurate results.
#' @param nmf_type Type of NMF to be used.
#'              Possible values:
#'              * 'regular'
#'              * 'variational_bayes'
#' The 'regular' method comes from the NMF package.
#' The 'variational_bayes' method comes from the ccfindR package.
#' This method uses bayesian inference, which makes it easier to determine the
#' mathematically optimal number of signatures.
#' @param single_core Boolean. If TRUE, it forces the NMF algorithm to
#' use only a single core. This can sometimes prevent issues.
#' Doesn't apply to variational-bayes NMF
#' @param fudge Small positive number that is used for the variational_bayes NMF.
#' Setting this to a small value like 0.0001 can prevent errors from occurring,
#' when extracting many signatures at once. In general, we recommend extracting
#' less signatures when errors occur, but this parameter can be used when that
#' is not an option.
#' Default = NULL.
#' @param seed Random seed used for the regular NMF, default = 123456
#'
#' @return Named list of mutation matrix, signatures and signature contribution
#'
#' @import NMF
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## This function is computationally intensive.
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' ## It's also possible to use a variational Bayes method.
#' ## It requires the ccfindR package to work.
#' # nmf_res <- extract_signatures(mut_mat, rank = 2, nmf_type = "variational_bayes")
#' @seealso
#' \code{\link{mut_matrix}}
#'
#' @export

extract_signatures <- function(mut_matrix, 
                               rank, 
                               nrun = 200, 
                               nmf_type = c("regular", "variational_bayes"), 
                               single_core = FALSE, 
                               fudge = NULL,
                               seed = 123456) {
  # Match argument
  nmf_type <- match.arg(nmf_type)

  # Add a small pseudocount to avoid features with zero counts.
  mut_matrix <- as.matrix(mut_matrix) + 0.0001

  # Make sure the rank_range is valid.
  if (!(rank > 0 & rank == round(rank))) {
    stop("Rank should be a positive integer", call. = FALSE)
  }

  if (ncol(mut_matrix) < max(rank)) {
    stop(paste0(
      "The rank should be smaller than the number of ",
      "samples in the input matrix."
    ), call. = FALSE)
  }

  if (nmf_type == "regular") {
    # Calculate NMF
    if (single_core){
      res <- NMF::nmf(mut_matrix, rank = rank, method = "brunet", nrun = nrun, seed = seed, .opt = "v-p")
    } else{
      res <- NMF::nmf(mut_matrix, rank = rank, method = "brunet", nrun = nrun, seed = seed)
    }
    # Find signatures and contribution of signatures
    signatures <- NMF::basis(res)
    contribution <- NMF::coef(res)
  } else {
    if (!requireNamespace("ccfindR", quietly = TRUE)) {
      stop(paste0(
        "Package 'ccfindR' is needed for variational_bayes to work. ",
        "Please either install it or use the regular NMF."
      ), call. = FALSE)
    }
    sc <- ccfindR::scNMFSet(count = mut_matrix)
    res <- ccfindR::vb_factorize(sc, ranks = rank, nrun = nrun, progress.bar = FALSE, verbose = 0, fudge = fudge)
    # estimate = ccfindR::vb_factorize(sc, ranks = 2:7, nrun = nrun, progress.bar = FALSE, verbose = 0)
    # plot(estimate)
    # optimal_rank(sb)
    signatures <- ccfindR::basis(res)[[1]]
    contribution <- ccfindR::coeff(res)[[1]]
  }
  # Reconstruct mutation matrix
  reconstructed <- signatures %*% contribution
  return(list(
    signatures = signatures,
    contribution = contribution,
    reconstructed = reconstructed
  ))
}

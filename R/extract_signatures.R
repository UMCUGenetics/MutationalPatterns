#' Extract mutational signatures from 96 mutation matrix using NMF
#'
#' Decomposes trinucleotide count matrix into signatures and contribution of
#' those signatures to the spectra of the samples/vcf files.
#'
#' @param mut_matrix 96 mutation count matrix 
#' @param rank Number of signatures to extract
#' @param nrun Number of iterations, default = 200
#' @return Named list of mutation matrix, signatures and signature contribution
#'
#' @import NMF
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## This function is computational intensive.
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' @seealso
#' \code{\link{mut_matrix}}
#'
#' @export

extract_signatures = function(mut_matrix, rank, nrun = 200)
{
    # Add a small pseudocount to avoid features with zero counts.
    mut_matrix = as.matrix(mut_matrix) + 0.0001

    # Make sure the rank_range is valid.
    if (!(rank > 0 & rank == round(rank)))
        stop("Rank should be a positive integer")

    if (ncol(mut_matrix) < max(rank))
        stop(paste( "The rank should be smaller than the number of",
                    "samples in the input matrix.") )

    # Calculate NMF
    res = nmf(mut_matrix, rank=rank, method="brunet", nrun=nrun, seed=123456)

    # Find signatures and contribution of signatures
    signatures = NMF::basis(res)
    contribution = NMF::coef(res)

    # Reconstruct mutation matrix
    reconstructed = signatures %*% contribution
    return(list(signatures = signatures,
                contribution = contribution,
                reconstructed = reconstructed))
}

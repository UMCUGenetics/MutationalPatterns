#' Find optimal nonnegative linear combination of mutation signatures to
#' reconstruct the mutation matrix.
#' 
#' Find the linear combination of mutation signatures that most closely
#' reconstructs the mutation matrix by solving the nonnegative least-squares
#' constraints problem.
#' 
#' @param mut_matrix 96 mutation count matrix (dimensions: 96 mutations
#' X n samples)
#' @param signatures Signature matrix (dimensions: 96 mutations
#' X n signatures)
#'
#' @return Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @importFrom pracma lsqnonneg
#'
#' @examples
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## You can download the signatures from the COSMIC website:
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
#'
#' ## Perform the fitting
#' fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
#'
#' @seealso
#' \code{\link{mut_matrix}}
#'
#' @export

fit_to_signatures = function(mut_matrix, signatures)
{
    # make sure dimensions of input matrix are correct
    if (dim(mut_matrix)[1] != dim(signatures)[1])
        stop(paste("Mutation matrix and signatures input have",
                   "different number of mutational features"))

    n_features = dim(mut_matrix)[1]
    n_samples = dim(mut_matrix)[2]
    n_signatures = dim(signatures)[2]
    lsq_contribution = matrix(NA, nrow=n_signatures, ncol=n_samples)
    lsq_reconstructed = matrix(NA, nrow=n_features, ncol=n_samples)

    # Process each sample
    for (i in 1:ncol(mut_matrix))
    {
        y = mut_matrix[,i]
        lsq = lsqnonneg(signatures, y)
        lsq_contribution[,i] = lsq$x
        lsq_reconstructed[,i] = signatures %*% as.matrix(lsq$x) 
    }

    # Add row and col names
    sample_names = colnames(mut_matrix)
    signature_names = colnames(signatures)
    mut_type_names = rownames(signatures)

    colnames(lsq_contribution) = sample_names
    rownames(lsq_contribution) = signature_names

    colnames(lsq_reconstructed) = sample_names
    rownames(lsq_reconstructed) = mut_type_names

    res = list(lsq_contribution, lsq_reconstructed)
    names(res) = c("contribution", "reconstructed")

    return(res)
}

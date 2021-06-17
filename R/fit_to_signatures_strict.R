#' Fit mutational signatures to a mutation matrix with less overfitting
#'
#' Refitting signatures with this function suffers less from overfitting.
#' The strictness of the refitting is dependent on 'max_delta'.
#' A downside of this method is that it might increase signature misattribution.
#' Different signatures might be attributed to similar samples.
#' You can use 'fit_to_signatures_bootstrapped()', to see if this is happening.
#' Using less signatures for the refitting will decrease this issue. Fitting
#' less strictly will also decrease this issue.
#'
#' Find a linear non-negative combination of mutation signatures that
#' reconstructs the mutation matrix. Signature selection (feature selection) is
#' done to reduce overfitting. This can be done via either a 'backwards'
#' (default) or 'best_subset' method. 
#' The 'backwards' method starts by achieving an optimal reconstruction via
#' 'fit_to_signatures'. The signature with the lowest contribution is then
#' removed and refitting is repeated. This is done in an iterative fashion. Each
#' time the cosine similarity between the original and reconstructed profile is
#' calculated.
#' The 'best_subset' method also starts by achieving an optimal reconstruction
#' via 'fit_to_signatures'. Signature refitting is then repeated for each
#' combination of n-1 signatures, where n is the number of signatures in the
#' signature matrix. The cosine similarity between the original and
#' reconstructed profile is calculated for each combination. The combination
#' with the highest cosine similarity is then chosen. This is done in an
#' iterative fashion for n-2, n-3, ect.
#' With both methods, iterations are stopped when the difference between two
#' iterations becomes more than 'max_delta'. The second-last set of signatures
#' is then used for a final refit.
#'
#' The 'best_subset' method can result in more accurate results than the 'backwards' method, however it becomes very slow when a large
#' amount of signatures are used for refitting. We recommend only using the 'best_subset' method when fitting a maximum of 10-15 signatures.
#' When using the 'best_subset' method a lower 'max_delta' should be used, as the expected differences in cosine similarity are reduced.
#'
#' @param mut_matrix Mutation count matrix (dimensions: x mutation types
#' X n samples)
#' @param signatures Signature matrix (dimensions: x mutation types
#' X n signatures)
#' @param max_delta The maximum difference in original vs reconstructed cosine similarity between two iterations.
#' @param method The method used to select signatures.
#' @return A list containing a fit_res object, similar to 'fit_to_signatures' and a list of ggplot graphs
#' that for each sample shows in what order the signatures were removed and how this affected the cosine similarity.
#'
#' @seealso \code{\link{mut_matrix}},
#' \code{\link{fit_to_signatures}},
#' \code{\link{fit_to_signatures_bootstrapped}}
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Get signatures
#' signatures <- get_known_signatures()
#'
#' ## Fit to signatures strict
#' strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
#'
#' ## fit_res similar to 'fit_to_signatures()'
#' fit_res <- strict_refit$fit_res
#'
#' ## list of ggplots that shows how the cosine similarity was reduced during the iterations
#' fig_l <- strict_refit$sim_decay_fig
#' 
#' 
#' ## Fit to signatures with the best_subset method
#' ## This can be more accurate than the standard backwards method, 
#' ## but can only be used with a limited amount of signatures.
#' ## Here we use only 5 signatures to reduce the runtime. 
#' ## In practice up to 10-15 signatures could be used.
#' best_subset_refit <- fit_to_signatures_strict(mut_mat, 
#'    signatures[,1:5], 
#'    max_delta = 0.002, 
#'    method = "best_subset"
#' )
#' 
fit_to_signatures_strict <- function(mut_matrix, signatures, max_delta = 0.004, method = c("backwards", "best_subset")) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  rowname <- . <- NULL
  
  # Match argument
  method <- match.arg(method)
  
  #Set colnames if absent, to prevent duplicate names later.
  if (is.null(colnames(mut_matrix))){
    colnames(mut_matrix) <- seq_len(ncol(mut_matrix))
  }

  # Remove signatures with zero contribution across samples
  fit_res <- fit_to_signatures(mut_matrix, signatures)
  sig_pres <- rowSums(fit_res$contribution) != 0
  my_signatures_total <- signatures[, sig_pres, drop = FALSE]

  # perform signature selection per sample
  all_results <- vector("list", ncol(mut_matrix))
  for (i in seq(1, ncol(mut_matrix))) {
    #my_signatures <- my_signatures_total
    mut_mat_sample <- mut_matrix[, i, drop = FALSE]

    if (method == "backwards"){
      results <- .strict_refit_backwards_selection_sample(mut_mat_sample, my_signatures_total, max_delta)
    } else if (method == "best_subset"){
      # Starts with all signatures
      results <- .strict_refit_best_subset_selection_sample(mut_mat_sample, signatures, max_delta)
    } else{
      stop("method needs to be either `backwards` or `best_subset`")
    }
    all_results[[i]] <- results
  }

  # Get decay figs and fit_res in separate lists
  decay_figs <- purrr::map(all_results, "sim_decay_fig")
  fit_res <- purrr::map(all_results, "fit_res")

  # Combine the contribution of all samples
  contribution <- purrr::map(fit_res, "contribution") %>%
    purrr::map(function(x) tibble::rownames_to_column(as.data.frame(x))) %>%
    purrr::reduce(dplyr::full_join, by = "rowname")

  # Fix signature order of contribution and add absent sigs to
  # keep the legend colors consistent for plotting.
  sig_ref <- tibble::tibble("rowname" = colnames(signatures))
  contribution <- dplyr::left_join(sig_ref, contribution, by ="rowname") %>% 
    as.data.frame()

  # Turn contribution into matrix and remove NAs
  rownames(contribution) <- contribution$rowname
  contribution <- contribution %>%
    dplyr::select(-rowname) %>%
    as.matrix()
  contribution[is.na(contribution)] <- 0

  # Combine the reconstructed of all samples
  reconstructed <- purrr::map(fit_res, "reconstructed") %>%
    do.call(cbind, .)

  # Combine all and return
  fit_res <- list("contribution" = contribution, "reconstructed" = reconstructed)
  results <- list("sim_decay_fig" = decay_figs, "fit_res" = fit_res)
  return(results)
}



#' Function to perform the strict signature refitting for a single sample with best subset selection
#'
#' @param mut_mat_sample mutation count matrix for a single sample
#' @param signatures signature matrix
#' @param max_delta The maximum difference in original vs reconstructed cosine similarity between two iterations.
#'
#' @return A list containing a fit_res object, similar to 'fit_to_signatures' and a ggplot graph
#' that for each sample shows in what order the signatures were removed and how this affected the cosine similarity.
#' 
#' @noRd
#'
.strict_refit_best_subset_selection_sample = function(mut_mat_sample, signatures, max_delta){
    
    # Determine the number of signatures.
    nsigs <- ncol(signatures)
    
    # Determine the starting similarity when all signatures are selected.
    old_res <- .strict_refit_best_subset_selection_sample_n(mut_mat_sample, signatures, nsigs)
    sims <- n_l <- vector("list", nsigs)
    sims[[1]] <- old_res$sim
    n_l[[1]] <- nsigs
    
    # For each n number of signatures select the best subset.
    for (n in seq(nsigs-1, 1)){
        
        # Find the best subset of n signatures
        new_res <- .strict_refit_best_subset_selection_sample_n(mut_mat_sample, signatures, n)
        sims[[nsigs-n+1]] <- new_res$sim
        n_l[[nsigs-n+1]] <- n
        
        # If the delta between two n's is below the cutoff switch to the next n.
        delta <- old_res$sim - new_res$sim
        if (delta <= max_delta){
            old_res <- new_res
        } else{
            break
        }
    }
    
    # Plot the decay in cosine similarity between decreasing subset sizes.
    sim_decay_fig <- .plot_sim_decay(sims, n_l, max_delta, "best_subset")
    
    
    # Do a final refit with the selected signatures.
    fit_res <- fit_to_signatures(mut_mat_sample, signatures[,old_res$best_subset, drop = FALSE])
    
    
    results <- list("sim_decay_fig" = sim_decay_fig, "fit_res" = fit_res)
    return(results)
}

#' Function to perform the strict signature refitting for a single sample with best subset selection for a single subset size
#'
#' @param mut_mat_sample Mutation count matrix for a single sample
#' @param signatures Signature matrix
#' @param n Number of signatures to use in the subset
#'
#' @return List containing the best subset of n signatures and the corresponding cosine similarity of 
#' the original vs reconstructed profile
#'
#' @noRd
#' 
.strict_refit_best_subset_selection_sample_n = function(mut_mat_sample, signatures, n){
    
    # Determine each possible combination of n signatures
    all_n_subsets <- combn(colnames(signatures), n, simplify = FALSE)
    
    # Perform a signature refit for each combination of n signatures.
    # And calculate the cosine similarity with the reconstructed profile.
    sims <- purrr::map_dbl(all_n_subsets, .strict_refit_best_subset_selection_single_set, mut_mat_sample, signatures)
    
    # Find which subest has the highest cosine similarity.
    best_subset <- all_n_subsets[[which.max(sims)]]
    best_sim <- max(sims)
    return(list("best_subset" = best_subset, "sim" = best_sim))
}


#' Function to determine the cosine similarity of the original vs reconstructed profile for a set of signatures.
#'
#' A refit is performed for a single sample, using a set of signatures.
#' The cosine similarity of the original vs reconstructed profile is then calculated.
#'
#' @param signature_names vector containing the names of the signatures to be used
#' @param mut_mat_sample Mutation count matrix for a single sample
#' @param signatures Signature matrix
#'
#' @return The cosine similarity of the original vs reconstructed profile
#'
#' @noRd
#' 
.strict_refit_best_subset_selection_single_set = function(signature_names, mut_mat_sample, signatures){
    
    # Perform a refit for a single signature set.
    signature_set <- signatures[,signature_names, drop = FALSE]
    fit_res <- fit_to_signatures(mut_mat_sample, signature_set)
    
    # Calculate the cosine similarity of the original vs the reconstructed profile.
    sim <- .get_cos_sim_ori_vs_rec(mut_mat_sample, fit_res)
    return(sim)
}


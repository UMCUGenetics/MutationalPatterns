#' Function to perform the strict signature refitting for a single sample with backwards selection
#'
#' @param mut_mat_sample mutation count matrix for a single sample
#' @param my_signatures signature matrix
#' @param max_delta The maximum difference in original vs reconstructed cosine similarity between two iterations.
#'
#' @return A list containing a fit_res object, similar to 'fit_to_signatures' and a ggplot graph
#' that for each sample shows in what order the signatures were removed and how this affected the cosine similarity.

#' @noRd
#'
.strict_refit_backwards_selection_sample = function(mut_mat_sample, my_signatures, max_delta){
    
    # Determine the number of signatures
    nsigs <- ncol(my_signatures)
    
    # Fit again
    fit_res <- fit_to_signatures(mut_mat_sample, my_signatures)
    sim <- .get_cos_sim_ori_vs_rec(mut_mat_sample, fit_res)
    
    # Keep track of the cosine similarity and which signatures are removed.
    sims <- vector("list", nsigs)
    sims[[1]] <- sim
    removed_sigs <- vector("list", nsigs)
    removed_sigs[[1]] <- "None"
    
    if (nsigs > 1){ # Only remove signatures if there is more than 1.
    # Sequentially remove the signature with the lowest contribution
        for (j in seq(2, nsigs)) {
            
            # Remove signature with the weakest relative contribution
            contri_order <- fit_res$contribution %>%
                prop.table(2) %>%
                rowSums() %>%
                order()
            weakest_sig_index <- contri_order[1]
            weakest_sig <- colnames(my_signatures)[weakest_sig_index]
            removed_sigs[[j]] <- weakest_sig
            signatures_sel <- my_signatures[, -weakest_sig_index, drop = FALSE]
            
            
            # Fit with new signature selection
            fit_res <- fit_to_signatures(mut_mat_sample, signatures_sel)
            sim_new <- .get_cos_sim_ori_vs_rec(mut_mat_sample, fit_res)
            
            if (is.nan(sim_new) == TRUE) {
                sim_new <- 0
                warning("New similarity between the original and the reconstructed 
                              spectra after the removal of a signature was NaN. 
                              It has been converted into a 0. 
                              This happened with the following fit_res:")
                print(fit_res)
            }
            sims[[j]] <- sim_new
            
            # Check if the loss in cosine similarity between the original vs reconstructed after removing the signature is below the cutoff.
            delta <- sim - sim_new
            if (delta <= max_delta) {
                my_signatures <- signatures_sel
                sim <- sim_new
            }
            else {
                break
            }
        }
    }
    
    # Plot how the cosine similarities decayed
    sim_decay_fig <- .plot_sim_decay(sims, removed_sigs, max_delta, "backwards")
    
    # Perform final fit on selected signatures
    fit_res <- fit_to_signatures(mut_mat_sample, my_signatures)
    
    # Add data of sample to list.
    results <- list("sim_decay_fig" = sim_decay_fig, "fit_res" = fit_res)
    return(results)
}

#' Get the cosine similarity between a reconstructed mutation matrix and the original
#'
#' @param mut_matrix mutation count matrix (dimensions: x mutation types
#' X n samples)
#' @param fit_res Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @return Cosine similarity
#' @noRd
#'
.get_cos_sim_ori_vs_rec <- function(mut_matrix, fit_res) {
    cos_sim_all <- cos_sim_matrix(mut_matrix, fit_res$reconstructed)
    cos_sim <- diag(cos_sim_all)
    mean_cos_sim <- mean(cos_sim)
    return(mean_cos_sim)
}


#' Plot decay in cosine similarity as signatures are removed.
#'
#' @param sims List of cosine similarities
#' @param removed_sigs List of iteratively removed signatures
#' @param max_delta The maximum difference in original vs reconstructed cosine similarity.
#' @param method The signature selection method that was used. Possible values:
#'              * 'backwards';
#'              * 'best_subset';
#' 
#' @import ggplot2
#' @importFrom magrittr %>%
#' @noRd
#' @return ggplot object
#'
.plot_sim_decay <- function(sims, removed_sigs, max_delta, method = c("backwards", "best_subset")) {
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    Removed_signatures <- Cosine_similarity <- NULL
    
    # Match argument
    method = match.arg(method)
    
    # Prepare data
    sims <- sims[!S4Vectors::isEmpty(sims)] %>%
        unlist()
    removed_sigs <- removed_sigs[!S4Vectors::isEmpty(removed_sigs)] %>%
        unlist()
    tb <- tibble::tibble(
        "Cosine_similarity" = sims,
        "Removed_signatures" = factor(removed_sigs, levels = removed_sigs)
    )
    
    # Determine if the final removed signature exceeded the cutoff.
    sims_l <- length(sims)
    col <- rep("low_delta", sims_l)
    if (sims_l > 1){ # Check if any signatures have been removed, before calculating the delta.
        final_delta <- sims[sims_l - 1] - sims[sims_l]
        if (final_delta > max_delta) {
            col[sims_l] <- "high_delta"
        }
    }
    # Set the x-axis label and theme
    if (method == "backwards"){
        xlab <- "Removed signatures"
        my_theme <- theme(
            axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 0.5),
            text = element_text(size = 12)
        )
    } else{
        xlab <- "Nr. signatures used"
        my_theme <- theme(text = element_text(size = 12))
    }
    
    # Create plot
    fig <- ggplot(data = tb, aes(x = Removed_signatures, y = Cosine_similarity, fill = col)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(
            limits = c("low_delta", "high_delta"),
            values = c("grey", "red"),
            guide = "none"
        ) +
        labs(
            x = xlab,
            y = paste0("Cosine similarity (max delta: ", max_delta, ")")
        ) +
        theme_classic() +
        my_theme
    return(fig)
}

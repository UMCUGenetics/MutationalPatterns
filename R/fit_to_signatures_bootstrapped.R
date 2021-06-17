#' Fit mutational signatures to a mutation matrix with bootstrapping
#'
#' @description
#' Bootstrapping the signature refitting shows how stable the refit is, when small changes are made to the
#' mutation matrix. You can be more confident in the refitting results, when the differences in signature
#' contributions are small between bootstrap iterations.
#'
#' @details
#' The mutation matrix is resampled 'n_boots' times.
#' Resampling is done per column (sample) with replacement.
#' The row weights are used as probabilities.
#' On each resampled matrix the 'fit_to_signatures()' or 'fit_to_signatures_strict()' function
#' is applied.
#' In the end a matrix is returned with the contributions for each bootstrap iteration.
#' Each row is a single bootstrap iteration from a single sample.
#' The method you choose determines how strict the signature refitting is.
#' The 'regular' and "regular_10+ methods often suffer from a lot of overfitting,
#' however this is less of an issue when you refit on an limited number of signatures.
#' The 'strict' method suffers less from overfitting, but can suffer from more
#' signature misattribution. The best method will depend on your data and
#' research question.
#'
#' @param mut_matrix mutation count matrix (dimensions: x mutation types
#' X n samples)
#' @param signatures Signature matrix (dimensions: x mutation types
#' X n signatures)
#' @param n_boots Number of bootstrap iterations.
#' @param max_delta The maximum difference in original vs reconstructed cosine similarity between two iterations.
#' Only used with method strict.
#' @param method The refitting method to be used.
#'               Possible values:
#'              * 'strict' Uses fit_to_signatures_strict with the default (backwards selection) method;
#'              * 'regular' Uses fit_to_signatures;
#'              * 'regular_10+' Uses fit_to_signatures, but removes signatures with less than 10 variants.;
#'              * 'strict_best_subset' Uses fit_to_signatures_strict with the 'best_subset' method;
#'              * 'strict_backwards' Uses fit_to_signatures_strict with the backwards selection method.
#'              This is the same as the 'strict' method;
#' @param verbose Boolean. If TRUE, the function will show how far along it is.
#'
#' @return A matrix showing the signature contributions across all the bootstrap iterations.
#' @export
#'
#' @seealso \code{\link{mut_matrix}},
#' \code{\link{fit_to_signatures_strict}},
#' \code{\link{fit_to_signatures_bootstrapped}}
#'
#' @importFrom magrittr %>%
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Get pre-defined signatures
#' signatures <- get_known_signatures()
#'
#' ## Fit to signatures with bootstrapping
#' ## Here we use a very low "n_boots" to reduce the runtime.
#' ## For real uses, a much higher value is required.
#' contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
#'   signatures,
#'   n_boots = 2,
#'   max_delta = 0.004
#' )
#'
#' ## Use the regular refit method
#' contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
#'   signatures,
#'   n_boots = 2,
#'   method = "regular"
#' )
fit_to_signatures_bootstrapped <- function(mut_matrix,
                                           signatures,
                                           n_boots = 1000,
                                           max_delta = 0.004,
                                           method = c("strict", "regular", "regular_10+", "strict_best_subset", "strict_backwards"),
                                           verbose = TRUE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  sigs <- . <- NULL

  # Match argument
  method <- match.arg(method)

  # Check enough mutations are present
  min_nr_muts <- colSums(mut_matrix) %>%
    min()
  if (min_nr_muts <= 10) {
    warning("At least one of your samples has less than 10 mutations. 
                This will negatively impact the signature refitting. 
                Please consider removing or pooling this sample.")
  }

  sig_names_tb <- tibble::tibble("sigs" = colnames(signatures))
  contri_list <- vector("list", n_boots)
  for (i in seq_len(n_boots)) {

    # Resample mut_mat
    mut_mat_resampled <- .resample_mut_mat(mut_matrix)

    # Perform refit method
    if (method == "strict" | method == "strict_backwards") {
      refit_out <- fit_to_signatures_strict(mut_mat_resampled, signatures, max_delta = max_delta)
      contri <- refit_out$fit_res$contribution
    }
    else if (method == "regular") {
      fit_res <- fit_to_signatures(mut_mat_resampled, signatures)
      contri <- fit_res$contribution
    }
    else if (method == "regular_10+") {
      fit_res <- fit_to_signatures(mut_mat_resampled, signatures)
      index <- rowSums(fit_res$contribution >= 10) != 0 # Check whether a signature has at least 10 mutations in a single sample
      contri <- fit_res$contribution[index, ]
    }
    else if (method == "strict_best_subset"){
      refit_out <- fit_to_signatures_strict(mut_mat_resampled, signatures, max_delta = max_delta, method = "best_subset")
      contri <- refit_out$fit_res$contribution
    }

    # Reformat contribution
    colnames(contri) <- paste(colnames(contri), i, sep = "_")
    contri_tb <- contri %>%
      as.data.frame() %>%
      tibble::rownames_to_column("sigs") %>%
      dplyr::left_join(sig_names_tb, ., by = "sigs") %>%
      dplyr::select(-sigs)

    # Add contribution to list
    contri_list[[i]] <- contri_tb

    if (i %% 10 == 0 & verbose == TRUE) {
      message(paste0("Performed ", i, " of ", n_boots, " iterations"))
    }
  }

  # Combine contribution from all samples
  contri_boots <- do.call(cbind, contri_list) %>%
    t()

  # Clean up format
  colnames(contri_boots) <- sig_names_tb$sigs
  contri_boots[is.na(contri_boots)] <- 0
  contri_boots <- contri_boots[, colSums(contri_boots) != 0, drop = FALSE]

  return(contri_boots)
}

#' Resample a mutation matrix
#'
#' The mutation matrix is resampled per column (sample).
#' Resampling is done with replacement using the row weights as propabilities.
#'
#' @param mut_matrix mutation count matrix (dimensions: x mutation types
#' X n samples)
#'
#' @return A resamples mutation matrix
#'
#' @noRd
#'
.resample_mut_mat <- function(mut_matrix) {
  mut_mat_resampled <- apply(mut_matrix, 2, function(x) {
    total_muts <- sum(x)
    sample_weights <- x / total_muts
    feature_rows <- sample(seq_along(x), total_muts, replace = TRUE, prob = sample_weights)
    row_counts <- table(feature_rows)
    index <- as.numeric(names(row_counts))
    x[index] <- as.vector(row_counts)
    x[-index] <- 0
    return(x)
  })
  return(mut_mat_resampled)
}

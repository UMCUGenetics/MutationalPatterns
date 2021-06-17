#' Count DBS contexts
#'
#' @details
#' Counts the number of DBS per COSMIC context from a GRanges or GRangesList object containing DBS variants.
#' This function applies the count_dbs_contexts_gr function to each gr in its input.
#' It then combines the results in a single tibble and returns this.
#'
#' @param vcf_list GRanges or GRangesList object containing DBS mutations in which the context was added with get_dbs_context.
#'
#' @return A tibble containing the number of DBS per COSMIC context per gr.
#'
#' @examples
#' ## Get a GRangesList or GRanges object with DBS contexts.
#' ## See 'dbs_get_context' for more info on how to do this.
#' grl_dbs_context <- readRDS(system.file("states/blood_grl_dbs_context.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' # Count the DBS contexts
#' count_dbs_contexts(grl_dbs_context)
#' @family DBS
#' @seealso \code{\link{get_dbs_context}}
#'
#' @export
count_dbs_contexts <- function(vcf_list) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  REF <- ALT <- NULL

  # Set possible ref and alt combis.
  categories <- DBS_CATEGORIES

  # Turn grl into list if needed.
  if (inherits(vcf_list, "CompressedGRangesList")) {
    vcf_list <- as.list(vcf_list)
  }

  # Count contexts per sample
  if (inherits(vcf_list, "list")) {
    counts_l <- purrr::map(vcf_list, .count_dbs_contexts_gr, categories)
    counts <- do.call(cbind, counts_l)
    colnames(counts) <- names(vcf_list)
  } else if (inherits(vcf_list, "GRanges")) {
    counts <- .count_dbs_contexts_gr(vcf_list, categories)
    colnames(counts) <- "My_sample"
  } else {
    .not_gr_or_grl(vcf_list)
  }
  counts <- cbind(categories, counts)
  counts[is.na(counts)] <- 0
  counts <- counts %>%
    tidyr::unite("muttype_total", REF, ALT) %>%
    tibble::column_to_rownames("muttype_total") %>%
    as.matrix()

  return(counts)
}



#' Count DBS contexts from a single GRanges object.
#'
#' @details
#' Counts the number of DBS per COSMIC context from a GRanges object containing DBS mutations.
#' The function is called by count_dbs_contexts
#'
#' @param gr GRanges object containing DBS mutations in which the context was added with 'get_dbs_context()'.
#' @param categories A tibble containing all possible DBS context categories
#'
#' @return A single column tibble containing the number of DBS per COSMIC context.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
.count_dbs_contexts_gr <- function(gr, categories) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  REF <- ALT <- NULL

  context <- cbind("REF" = as.vector(.get_ref(gr)), "ALT" = as.vector(unlist(.get_alt(gr))))
  counts <- context %>%
    tibble::as_tibble() %>%
    dplyr::group_by(REF, ALT) %>%
    dplyr::summarise(count = dplyr::n())
  
  if (sum(!counts$REF %in% categories$REF) > 0 | sum(!counts$ALT %in% categories$ALT)){
    stop(paste0("There are some REF or ALT bases, that do not belong to ", 
                "any of the categories. \n",
                "Did you forget to use 'get_dbs_context()'?"), call. = FALSE)
  }
  
  
  counts_full <- dplyr::left_join(categories, counts, by = c("REF", "ALT")) %>%
    dplyr::select(-REF, -ALT)
  return(counts_full)
}

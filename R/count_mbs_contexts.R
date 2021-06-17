#' Count MBS variants grouped by length.
#'
#' @details
#' Counts the number of mbs grouped by length from a GRanges or GRangesList object containing mbs variants.
#' This is used, since a COSMIC context has to our knowledge not yet been defined.
#' This function applies the count_mbs_contexts_gr function to each gr in its input.
#' It then combines the results in a single tibble and returns this.
#'
#' @param vcf_list GRanges or GRangesList object containing mbs variants.
#'
#' @return A tibble containing the number of MBS per MBS length per gr.
#'
#' @examples
#' ## Get a GRangesList or GRanges object with mbs variants.
#' mbs_grl <- readRDS(system.file("states/blood_grl_mbs.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' # Count the MBSs
#' count_mbs_contexts(mbs_grl)
#' @family MBS
#'
#' @export
count_mbs_contexts <- function(vcf_list) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  size <- NULL

  categories <- MBS_CATEGORIES

  # Turn grl into list if needed.
  if (inherits(vcf_list, "CompressedGRangesList")) {
    vcf_list <- as.list(vcf_list)
  }

  # Count contexts per sample
  if (inherits(vcf_list, "list")) {
    counts_l <- purrr::map(vcf_list, .count_mbs_contexts_gr, categories)
    counts <- do.call(cbind, counts_l)
    colnames(counts) <- names(vcf_list)
  } else if (inherits(vcf_list, "GRanges")) {
    counts <- .count_mbs_contexts_gr(vcf_list, categories)
    colnames(counts) <- "My_sample"
  } else {
    .not_gr_or_grl(vcf_list)
  }
  counts <- cbind(categories, counts)

  # Turn output into matrix
  counts <- counts %>%
    tibble::column_to_rownames("size") %>%
    as.matrix()
  return(counts)
}

#' Count MBS grouped by length from a single GRanges object.
#'
#' @details
#' Counts the number of MBS per MBS length from a GRanges object containing mbs variants.
#' The function is called by count_mbs_contexts
#'
#' @param gr GRanges object containing mbs variants.
#' @param categories A tibble containing all possible mbs size categories
#'
#' @return A single column tibble containing the number of MBS per MBS length
#'
#' @importFrom magrittr %>%
#' @noRd
.count_mbs_contexts_gr <- function(gr, categories) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  count <- size <- . <- NULL

  # Create count table
  counts_tb <- gr %>% # Determine different sizes
    .get_ref() %>%
    BiocGenerics::width() %>%
    tibble::enframe(value = "size") %>%
    dplyr::select(-name) %>%
    dplyr::mutate(
      size = ifelse(size >= 10, "10+", size),
      size = as.character(size)
    ) %>%
    dplyr::group_by(size) %>% # Count number muts per size
    dplyr::summarise(count = dplyr::n(), .groups = "drop_last") %>%
    dplyr::left_join(categories, ., by = "size") %>% # Add to possible categories
    dplyr::mutate(count = ifelse(is.na(count), 0, count)) %>%
    dplyr::select(-size)

  return(counts_tb)
}

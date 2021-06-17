
#' Count indel contexts
#'
#' @details
#' Counts the number of indels per COSMIC context from a GRanges or GRangesList object containing indel mutations.
#' This function applies the count_indel_contexts_gr function to each gr in its input.
#' It then combines the results in a single tibble and returns this.
#'
#' @param vcf_list GRanges or GRangesList object containing indel mutations in which the context was added with get_indel_context.
#'
#' @return A tibble containing the number of indels per COSMIC context per gr.
#'
#' @examples
#' ## Get a GRangesList or GRanges object with indel contexts.
#' ## See 'indel_get_context' for more info on how to do this.
#' grl_indel_context <- readRDS(system.file("states/blood_grl_indel_context.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' # Count the indel contexts
#' count_indel_contexts(grl_indel_context)
#' @family Indels
#'
#' @seealso \code{\link{get_indel_context}}
#'
#' @export
count_indel_contexts <- function(vcf_list) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  muttype <- muttype_sub <- NULL

  categories <- INDEL_CATEGORIES

  # Turn grl into list if needed.
  if (inherits(vcf_list, "CompressedGRangesList")) {
    vcf_list <- as.list(vcf_list)
  }

  # Count contexts per sample
  if (inherits(vcf_list, "list")) {
    counts_l <- purrr::map(vcf_list, .count_indel_contexts_gr, categories)
    counts <- do.call(cbind, counts_l)
    colnames(counts) <- names(vcf_list)
  } else if (inherits(vcf_list, "GRanges")) {
    counts <- .count_indel_contexts_gr(vcf_list, categories)
    colnames(counts) <- "My_sample"
  } else {
    .not_gr_or_grl(vcf_list)
  }
  counts <- cbind(categories, counts)
  counts[is.na(counts)] <- 0
  counts <- counts %>%
    tidyr::unite("muttype_total", muttype, muttype_sub) %>%
    tibble::column_to_rownames("muttype_total") %>%
    as.matrix()

  # counts = dplyr::as_tibble(counts)
  # counts$muttype = factor(counts$muttype, levels = unique(counts$muttype))
  return(counts)
}

#' Count indel contexts from a single GRanges object.
#'
#' @details
#' Counts the number of indels per COSMIC context from a GRanges object containing indel mutations.
#' The function is called by count_indel_contexts
#'
#' @param gr GRanges object containing indel mutations in which the context was added with get_indel_context.
#' @param categories A tibble containing all possible indel context categories
#'
#' @return A single column tibble containing the number of indels per COSMIC context.
#'
#' @importFrom magrittr %>%
#' @noRd
#'
.count_indel_contexts_gr <- function(gr, categories) {
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  muttype <- muttype_sub <- NULL

  # Check gr is not empty
  if (length(gr) == 0) {
    categories <- categories %>%
      dplyr::mutate(count = 0) %>%
      dplyr::select(-muttype, -muttype_sub)
    return(categories)
  }

  # Check context has previously been set.
  gr_colnames <- colnames(mcols(gr))
  if (!all(c("muttype", "muttype_sub") %in% gr_colnames)) {
    stop("The GRanges object does not contain the columns `muttype`` and `muttype_sub`.
             Did you forget to run `get_indel_context`?", call. = FALSE)
  }


  # Classify the number of repeat units/ homopolymer length / microhomology length
  # to either 5+ or 6+ depending on whether the indel is a ins or del.
  id_context <- dplyr::tibble("muttype" = gr$muttype, "muttype_sub" = gr$muttype_sub) %>%
    dplyr::mutate(
      muttype_sub = ifelse(muttype_sub >= 6, "6+", muttype_sub),
      muttype_sub = ifelse(grepl("insertion|microhomology", muttype) & muttype_sub >= 5,
        "5+", muttype_sub
      ),
      muttype_sub = as.character(muttype_sub)
    ) # Ensures column type for later joining

  # Classify large indels as size 5+
  id_context$muttype = .set_large_indels_as_5plus(id_context$muttype, gr)
  
  id_context_count <- id_context %>%
    dplyr::group_by(muttype, muttype_sub) %>%
    dplyr::summarise(count = dplyr::n())
  id_context_count_full <- dplyr::left_join(categories,
    id_context_count,
    by = c("muttype", "muttype_sub")
  ) %>%
    dplyr::select(-muttype, -muttype_sub)
  # colnames(id_context_count_full)[3] = name
  return(id_context_count_full)
}


#' Classify large indels as size 5+
#' Indels with a size that is more or equal to 5,
#' are set to size 5+
#'
#' @param muttype A vector containing the contexts of the mutations
#' @param gr GRanges object containing indel mutations in which the context was added with get_indel_context.
#'
#' @return A modified version of the vector containing the contexts of the mutations
#' 
#' @importFrom magrittr %>%
#' @noRd
.set_large_indels_as_5plus = function(muttype, gr){
  
  # Determine mutation sizes
  ref_sizes <- gr %>%
    .get_ref() %>%
    width()
  alt_sizes <- gr %>%
    .get_alt() %>%
    unlist() %>%
    width()
  mut_size <- abs(alt_sizes - ref_sizes)
  
  # Change classification of large mutations into 5+
  mut_size_f <- mut_size >= 5
  muttype <- ifelse(mut_size_f,
                               gsub("[0-9]+bp",
                                    "5+bp",
                                    muttype,
                                    perl = TRUE
                               ),
                               muttype
  )
  return(muttype)
}
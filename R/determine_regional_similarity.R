#' Determine regional mutation pattern similarity
#'
#' Calculate the cosine similarities between the global mutation profile and the
#' mutation profile of smaller genomic windows, using a sliding window approach.
#' Regions with a very different mutation profile can be identified in this way.
#' This function generally requires many mutations (~100,000) to work properly.
#'
#' @details First a global mutation matrix is calculated using all mutations.
#'   Next, a sliding window is used. This means that we create a window
#'   containing the first x mutations. The cosine similarity, between the
#'   mutation profiles of this window and the global mutation matrix, is then
#'   calculated. The window then slides y mutations to the right and the cosine
#'   similarity is again calculated. This process is repeated until the final
#'   mutation on a chromosome is reached. This process is performed separately
#'   per chromosome. Windows that span a too large region of the genome are
#'   removed, because they are unlikely to contain biologically relevant
#'   information.
#'   
#'   The number of mutations that the window slides to the right in each step is
#'   called the stepsize. The best stepsize depends on the window size. In
#'   general, we recommend setting the stepsize between 25% and 100% of the
#'   window size.
#'
#'   The analysis can be performed for trinucleotides contexts, for a larger
#'   context, or for just the base substitutions. A smaller context might miss
#'   detailed differences in mutation profiles, but is also less noisy. We
#'   recommend using a smaller extension when analyzing small datasets.
#'
#'   It's possible to correct for the oligonucleotide frequency of the windows.
#'   This is done by calculating the cosine similarity of the oligonucleotide
#'   frequency between each window and the genome. The cosine similarity of the
#'   mutation profiles is then divided by the oligonucleotide similarity. This
#'   ensures that regions with an abnormal oligonucleotide frequency don't show
#'   up as having a very different profile. The oligonucleotide frequency
#'   correction slows down the function, so we advise the user to keep it off
#'   for exploratory analyses and to only turn it on to validate interesting
#'   results.
#'
#'   By default the mutations in a window are subtracted from the global
#'   mutation matrix, before calculating the cosine similarity. This increases
#'   sensitivity, but could also decrease specificity. This subtraction can be
#'   turned of with the 'exclude_self_mut_mat' argument.
#'
#' @param vcf GRanges object
#' @param ref_genome BSgenome reference genome object
#' @param chromosomes Vector of chromosome/contig names of the reference genome
#'   to be plotted.
#' @param window_size The number of mutations in a window. (Default: 100)
#' @param stepsize The number of mutations that a window slides in each step.
#'   (Default: 25)
#' @param extension The number of bases, that's extracted upstream and
#'   downstream of the base substitutions, to create the mutation matrices.
#'   (Default: 1).
#' @param oligo_correction Boolean describing whether oligonucleotide frequency
#'   correction should be applied. (Default: FALSE)
#' @param exclude_self_mut_mat Boolean describing whether the mutations in a
#'   window should be subtracted from the global mutation matrix. (Default:
#'   TRUE)
#' @param max_window_size_gen The maximum size of a window before it is removed.
#'   (Default: 20,000,000)
#' @param verbose Boolean determining the verbosity of the function. (Default: FALSE)
#'
#' @return A "region_cossim" object containing both the cosine similarities and
#'   the settings used in this analysis.
#'
#' @seealso \code{\link{plot_regional_similarity}}
#' @family regional_similarity
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#'
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## We pool all the variants together, because the function doesn't work well
#' ## with a limited number of mutations. Still, in practice we recommend to use
#' ## more mutations that in this example.
#' gr = unlist(grl)
#'
#' ## Specifiy the chromosomes of interest.
#' chromosomes <- names(genome(gr)[1:3])
#'
#'  ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Determine the regional similarities. Here we use a small window size to make the function work.
#' ## In practice, we recommend a larger window size.
#' regional_sims = determine_regional_similarity(gr,
#'   ref_genome, 
#'   chromosomes,
#'   window_size = 40,
#'   stepsize = 10,
#'   max_window_size_gen = 40000000
#' )
#' 
#' ## Here we use an extensiof of 0 to reduce noise.
#' ## We also turned verbosity on, so you can see at what step the function is.
#' ## This can be useful on large datasets.
#' regional_sims_0_extension = determine_regional_similarity(gr,
#'   ref_genome, 
#'   chromosomes,
#'   window_size = 40,
#'   stepsize = 10,
#'   extension = 0,
#'   max_window_size_gen = 40000000,
#'   verbose = TRUE
#' )
determine_regional_similarity <- function(vcf,
                                          ref_genome,
                                          chromosomes,
                                          window_size = 100,
                                          stepsize = 25,
                                          extension = 1,
                                          oligo_correction = FALSE,
                                          exclude_self_mut_mat = TRUE,
                                          max_window_size_gen = 20000000,
                                          verbose = FALSE){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    chr <- . <- window_pos <- NULL
    
    # Get the reference genome
    tryCatch(
        error = function(cnd) {
            stop("Please provide the name of a BSgenome object.", call. = FALSE)
        },
        {
            ref_genome <- BSgenome::getBSgenome(ref_genome)
        }
    )

    # Check the input is a GRanges object.
    if (!inherits(vcf, "GRanges")) {
       .not_gr(vcf)
    }
    # Determine global mutation matrix.
    # This is done before preprocessing the gr, 
    # so that subsetting a single chromosome won't affect the global mutation matrix.
    mut_mat_total <- mut_matrix(vcf, ref_genome, extension)
    
    if (verbose){
        message("Created main mutation matrix")
    }
    
    # Preprocess the gr
    GenomeInfoDb::seqlevels(vcf, pruning.mode = "coarse") <- chromosomes
    vcf <- BiocGenerics::sort(vcf)
    
    if (verbose){
        message("Sorted the gr")
    }
    
    # get chromosome lengths of reference genome
    chr_lengths <- GenomeInfoDb::seqlengths(vcf)
    
    # Check for missing chromosome lengths
    if (sum(is.na(GenomeInfoDb::seqlengths(vcf))) > 1) {
        stop(paste(
            "Chromosome lengths missing from vcf object.\n",
            "Likely cause: contig lengths missing from the header of your vcf file(s).\n",
            "Please evaluate: seqinfo(vcf)\n",
            "To add seqlengths to your vcf GRanges object use: seqlengths(vcf) <- my_chr_lengths" 
        ), call. = FALSE)
    }
    
    
    # Split the vcf per chromosome and remove unused chromosome levels
    grl <- split(vcf, seqnames(vcf), drop = TRUE)
    GenomeInfoDb::seqlevels(grl) <- GenomeInfoDb::seqlevelsInUse(grl)
    muts_per_chr <- S4Vectors::elementNROWS(grl)
    
    # Determine the global nucleotide context if requested.
    if (oligo_correction == T){
        all_seq <- BSgenome::getSeq(ref_genome, chromosomes)
        global_oligo_context <- .get_oligo_contexts(all_seq, chromosomes, extension) %>% 
            rowSums() %>%
            as.matrix()
    } else{
        global_oligo_context <- NA
    }
    # Calculate the cosine similarities between local windows and all global mutations per chromosome.
    gr_l <- as.list(grl)
    sim_tb_l <- purrr::map(gr_l, ~.determine_regional_similarity_chr(.x, 
                                                                     ref_genome = ref_genome,
                                                                     window_size = window_size,
                                                                     stepsize,
                                                                     extension,
                                                                     mut_mat_total,
                                                                     global_oligo_context,
                                                                     oligo_correction,
                                                                     exclude_self_mut_mat,
                                                                     max_window_size_gen,
                                                                     verbose))
    
    # Extract the cosine similaries
    sim_tb <- sim_tb_l %>%
        purrr::map("sim_tb") %>%
        do.call(rbind, .) %>%
        dplyr::mutate(chr = factor(chr, levels = unique(chr)), window_pos_mb = window_pos / 1000000)
    mean_window_size <- mean(sim_tb$window_sizes)
    
    # Extract the mutation positions.
    pos_tb <- sim_tb_l %>%
        purrr::map("pos_tb") %>%
        do.call(rbind, .) %>%
        dplyr::mutate(chr = factor(chr, levels = unique(chr)), pos_mb = pos / 1000000)
    
    # Combine all results and settings in a single object.
    sims <- new("region_cossim",
               "sim_tb" = sim_tb,
               "pos_tb" = pos_tb,
               "chr_lengths" = chr_lengths,
               "window_size" = window_size,
               "max_window_size_gen" = max_window_size_gen,
               "ref_genome" = ref_genome,
               "muts_per_chr" = muts_per_chr,
               "mean_window_size" = mean_window_size,
               "stepsize" = stepsize,
               "extension" = extension,
               "chromosomes" = chromosomes,
               "exclude_self_mut_mat" = exclude_self_mut_mat)
    return(sims)
}

#' Calculate the number of mutations per oligonucleotide context per window
#'
#' @param seqs DNAStringSet of windows.
#' @param window_names A vector describing the names of each window (Default: my_window).
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions, to create the mutation matrices. 
#' 
#' @return A dataframe showing the number of mutations per oligonucleotide context per window.
#' @noRd
#' @importFrom magrittr %>%
#'
.get_oligo_contexts <- function(seqs, window_names = "my_window", extension){
    
    
    # Determine oligo nucleotide frequencies.
    oligo_width = 1 + (2 * extension)
    oligo_contexts = Biostrings::oligonucleotideFrequency(seqs, width = oligo_width)
    
    # Get the results in a tibble.
    if (inherits(oligo_contexts, "integer")){
        oligo_contexts_t <- tibble::enframe(oligo_contexts, name = "type", value = window_names)
    } else if (inherits(oligo_contexts, "matrix")){
        oligo_contexts_t <- oligo_contexts %>% 
            t() %>% 
            as.data.frame() %>% 
            tibble::rownames_to_column("type")
    }
    
    # Reverse complement the G and A mutations
    types <- oligo_contexts_t$type
    substitution_i = ceiling(oligo_width/2)
    bases <- types %>% 
        stringr::str_sub(start = substitution_i, end = substitution_i)
    context_to_reverse_i <- bases %in% c("G", "A")
    types_reversed <- types[context_to_reverse_i] %>% 
        Biostrings::DNAStringSet() %>% 
        Biostrings::reverseComplement() %>%
        as.character()
    oligo_contexts_t$type[context_to_reverse_i] <- types_reversed
    
    # Sum the contexts together
    oligo_contexts_t <- oligo_contexts_t %>% 
        dplyr::group_by(type) %>% 
        dplyr::summarise_all(sum) %>% 
        tibble::column_to_rownames("type")
    return(oligo_contexts_t)
}

#' Determine regional mutation pattern similarity for one chromosome
#'
#' @param gr GRanges object
#' @param ref_genome BSgenome reference genome object
#' @param window_size The number of mutations in a window.
#' @param stepsize The number of mutations that a window slides in each step.
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions, to create the mutation matrices. 
#' @param mut_mat_total Matrix containing mutation counts for all mutations.
#' @param global_oligo_context Matrix containing the global oligonucleotide frequencies. This is only used when,
#' correcting for oligonucleotide frequency.
#' @param oligo_correction Boolean describing whether oligonucleotide frequency correction should be applied.
#' @param exclude_self_mut_mat Boolean describing whether the mutations in a window should be subtracted from the global mutation matrix.
#' @param max_window_size_gen The maximum size of a window before it is removed.
#' @param verbose Boolean determining the verbosity of the function. (Default: FALSE)
#'
#' @return A list containing both the calculated similarities of the windows and the mutation positions.
#' @noRd
#' @importFrom magrittr %>%
#'
.determine_regional_similarity_chr <- function(gr, 
                                 ref_genome,
                                 window_size, 
                                 stepsize,
                                 extension,
                                 mut_mat_total, 
                                 global_oligo_context, 
                                 oligo_correction,
                                 exclude_self_mut_mat,
                                 max_window_size_gen,
                                 verbose){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    cossim <- . <- window_context_sim <- NULL
    
    # Determine chromosome
    chr <- GenomeInfoDb::seqlevelsInUse(gr)
   
    # Get base positions
    pos <- start(gr)
    pos_tb <- tibble::tibble("chr" = chr, "pos" = pos)
    
    # Check if window size is not too large
    if (window_size > length(gr)){
        message(paste0("Window size is larger than the number of mutations for: ", chr, ". Returning empty for this chromosome."))
        sim_tb <- dplyr::as_tibble(data.frame(matrix(nrow=0,ncol = 7)))
        colnames(sim_tb) <- c("cossim", "chr", "window_pos", "window_begin", "window_end", "window_sizes", "window_sizes_mb")
        return(list("sim_tb" =  sim_tb, "pos_tb" = pos_tb))
    }
    
    # Calculate window sizes in bp
    window_size_gen <- dplyr::lead(pos, n = window_size-1) - pos + 1
    window_size_gen <- window_size_gen[!is.na(window_size_gen)]
    
    # Determine starting index of each window
    window_i <- seq(1, length(window_size_gen), by = stepsize)
    window_size_gen <- window_size_gen[window_i]
    
    # Remove windows that are larger than the cutoff
    window_size_gen_f <- window_size_gen <= max_window_size_gen
    window_i <- window_i[window_size_gen_f]
    window_sizes <- window_size_gen[window_size_gen_f]
    
    # Return early if no windows were found
    if (length(window_i) == 0){
        message(paste0("No windows following the requirements were found in chr: ", chr))
        sim_tb <- dplyr::as_tibble(data.frame(matrix(nrow=0,ncol = 7)))
        colnames(sim_tb) <- c("cossim", "chr", "window_pos", "window_begin", "window_end", "window_sizes", "window_sizes_mb")
        return(list("sim_tb" =  sim_tb, "pos_tb" = pos_tb))
    }
    
    # Determine middle of windows
    if (window_size %% 2 == 0){
        pos_i <- window_i + window_size / 2 - 1
    } else{
        pos_i <- window_i + window_size/2-0.5
    }
    
    window_pos <- pos[pos_i]
    
    if (verbose){
        message(paste0("Determined window locations for chromosome: ", chr))
    }
    
    # Create an index containing the position of all mutations of window 1, followed by window 2, ect.
    nr_windows <- length(window_i)
    window_i_end <- window_i + window_size - 1
    seq_l <- mapply(seq, from = window_i, to = window_i_end, SIMPLIFY = FALSE)
    all_seq <- do.call(c, seq_l)
    
    
    # Determine the context of each mutation. 
    # Then use the 'all_seq' position index to repeat the contexts for muts that occur in more than 1 window.
    # This is then used to create the mut matrix.
    types <- type_context(gr, ref_genome, extension)
    types_allwindows <- list("types" = types$types[all_seq], "context" = types$context[all_seq])
    mut_mat <- mut_96_occurrences(types_allwindows, rep(window_size, nr_windows))
    
    if (verbose){
        message(paste0("Finished the mutation matrix for chromosome: ", chr))
    }
    
    # Determine cosine similarity of local mut mats with the global.
    if (exclude_self_mut_mat){
        cos_sims <- .cos_sim_matrix_exclself(mut_mat, mut_mat_total)  
    } else{
        cos_sims <- cos_sim_matrix(mut_mat, mut_mat_total)
    }
    colnames(cos_sims) <- "cossim"
    
    # Calculate some final locations and return tibbles.
    window_begin <- pos[window_i]
    window_end <- pos[window_i_end]
    
    # Determine trinucleotide context of windows. And its similarity to the entire genome.
    if (oligo_correction){
        
        
        # This is done in a loop to prevent vector memory exhaustion.
        context_sim <- purrr::map(seq(1, nr_windows, by = 10000), function(i){
            end <- min(c(i + 9999, nr_windows))
            
            # Get the sequence of each window.
            window_seqs <- Biostrings::getSeq(ref_genome, rep(chr, length(window_begin[i:end])), window_begin[i:end], window_end[i:end])
            
            # Get the tri context
            context_windows <- .get_oligo_contexts(window_seqs, extension = extension)
            
            # Determine the cosine similarity in nucleotide contexts of the windows with the global context.
            context_sim <- cos_sim_matrix(context_windows, global_oligo_context)
        }) %>% do.call(rbind, .)

        
        sim_tb <- tibble::tibble("cossim" = cos_sims[,"cossim"],
                        "window_context_sim" = context_sim[,1], chr, window_pos, window_begin, window_end, window_sizes) %>% 
            dplyr::mutate(window_sizes_mb = window_sizes / 1000000, "corrected_cossim" = cossim / window_context_sim)
    } else{
        sim_tb <- tibble::tibble("cossim" = cos_sims[,"cossim"], chr, window_pos, window_begin, window_end, window_sizes) %>% 
            dplyr::mutate(window_sizes_mb = window_sizes / 1000000)
    }
    
    if (verbose){
        message(paste0("Finished calculating for chromosome: ", chr))
    }
    return(list("sim_tb" =  sim_tb, "pos_tb" = pos_tb))
}


#' Calculate cosine similarity, excluding mutations from window.
#' 
#' For each window the mutations from that window are subtracted from the global
#' mutation matrix. Then the cosine similarity between the window and the
#' subtracted global matrix are calculated.
#'
#' @param mut_mat Matrix containing mutation counts.
#' @param mut_mat_total Matrix containing mutation counts for all mutations.
#'
#' @return Matrix containing cosine similarities.
#' @noRd
#'
.cos_sim_matrix_exclself <- function(mut_mat, mut_mat_total){
    nr_windows <- ncol(mut_mat)
    
    # Loop over all windows
    cossims <- vector("list", nr_windows)
    for (i in seq_len(nr_windows)){
        window_mut_mat <- mut_mat[,i]
        
        # Remove window from global mutation matrix
        mut_mat_total_exclself <- mut_mat_total[,1] - window_mut_mat
        
        # Calculate cosine similarities
        cossim <- cos_sim(window_mut_mat, mut_mat_total_exclself)
        cossims[[i]] <- cossim
    }
    
    # Combine all results
    cossims <- do.call(rbind, cossims)
    return(cossims)
}

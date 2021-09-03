#' Plot regional similarity
#' 
#' Plot the cosine similarity of the mutation profiles of small genomic windows
#' with the rest of the genome.
#' 
#' @details
#' Each dot shows the cosine similarity between the mutation profiles of a
#' single window and the rest of the genome. A region with a different mutation
#' profile will have a lower cosine similarity. The dots are colored based on
#' the sizes in mega bases of the windows. This size is the distance between the
#' first and last mutations in a window. The locations of the mutations can be
#' plotted on the bottom of the figure. The cosine similarity can be plotted
#' both with and without oligonucleotide frequency correction. This can be done
#' for all chromosomes at once or separate plots can be made per chromosome.
#'
#' @param region_cossim A region_cossim object.
#' @param per_chrom Boolean. Determines whether to create a separate plot per chromosome. (Default: FALSE)
#' @param oligo_correction Boolean describing whether the oligonucleotide
#'   frequency corrected cosine similarities should be plotted. If no correction
#'   has been applied then the regular cosine similarities will be plotted.
#'   (Default: TRUE)
#' @param max_cossim Maximum cosine similarity for a window to be considered
#'   an outlier. Any window with a lower cosine similarity is given a different
#'   color. (Default: NA)
#' @param title Optional plot title. (Default: NA). When the default option is
#'   used, the number of mutations per window and the step size are shown.
#' @param plot_rug Add a bottom rug to the plot, depicting the location of
#'   the mutations. (Default: FALSE)
#' @param x_axis_breaks Vector of custom x-axis breaks. (Default: NA)
#'
#' @return ggplot2 object
#' 
#' @export
#' @seealso
#' \code{\link{determine_regional_similarity}}
#' @family regional_similarity
#' 
#' @importFrom magrittr %>%
#' 
#' @examples
#' 
#' ## See the 'determine_regional_similarity()' example for how we obtained the
#' ## following data:
#' regional_sims <- readRDS(system.file("states/regional_sims.rds",
#'   package = "MutationalPatterns"
#' ))
#' 
#' ## Plot the regional similarity
#' plot_regional_similarity(regional_sims)
#' 
#' ## Plot outlier samples with a different color.
#' ## The value of 0.5 that is used here is arbitrarily chosen
#' ## and should in practice be based on the data.
#' plot_regional_similarity(regional_sims, max_cossim = 0.5)
#' 
#' ## Plot samples per chromosome
#' fig_l = plot_regional_similarity(regional_sims, per_chrom = TRUE)
#' 
#' ## Plot without a title
#' plot_regional_similarity(regional_sims, title = "")
#' 
#' ## Add a rug to the plot, that shows the location of the mutations.
#' plot_regional_similarity(regional_sims, plot_rug = FALSE)
#' 
#' ## Use custom x axis breaks
#' plot_regional_similarity(regional_sims, x_axis_breaks = c(50, 150))
#' 
plot_regional_similarity <- function(region_cossim, 
                                    per_chrom = FALSE, 
                                    oligo_correction = TRUE, 
                                    max_cossim = NA,
                                    title = NA,
                                    plot_rug = FALSE,
                                    x_axis_breaks = NA){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    chr <- sims <- metadata <- NULL
    
    #Retrieve chromosome lengths.
    max_chr_length <- max(region_cossim@chr_lengths / 1000000)
    
    #Set x-axis breaks
    if (.is_na(x_axis_breaks)){
        x_axis_breaks <- .get_x_axis_breaks(region_cossim, per_chrom)
    }
    
    #Get pos_tb and sim_tb.
    pos_tb <- region_cossim@pos_tb
    sim_tb <- get_sim_tb(region_cossim)
    
    # Set oligo_correction to FALSE if it has not been performed
    if (!("corrected_cossim" %in% colnames(sim_tb))){
        oligo_correction = FALSE
    }
    
    #Remove chr from the chrom names, so they are shorter
    levels(pos_tb$chr) <- gsub("chr|chromosome_|chromosome|group|group_|chrom", "", levels(pos_tb$chr))
    levels(sim_tb$chr) <- gsub("chr|chromosome_|chromosome|group|group_|chrom", "", levels(sim_tb$chr))
    sim_tb <- sim_tb %>% 
        dplyr::mutate(chr = factor(chr, levels = levels(pos_tb$chr)))
    

    #Set point size
    windows_per_bp <- nrow(sim_tb) / sum(region_cossim@chr_lengths)
    point_size <- 0.0000008 / windows_per_bp
    
    if (per_chrom == TRUE){
        point_size <- point_size * 5
    }
    
    if (point_size > 1.5){
        point_size <- 1.5
    } else if (point_size < 0.03){
        point_size <- 0.03
    }
    
    chrom_ends <- .calc_chrom_ends(region_cossim)
    
    #Create plots
    if (per_chrom == FALSE){
        fig <- .plot_regional_similarity_gg(sim_tb, pos_tb, chrom_ends, x_axis_breaks, point_size, region_cossim@window_size,
                                    region_cossim@stepsize, oligo_correction, max_cossim = max_cossim, title, plot_rug)
        return(fig)   
    } else{
        sim_tb_l <- split(sim_tb, sim_tb$chr)
        pos_tb_l <- split(pos_tb, pos_tb$chr)
        chrom_ends_l <- split(chrom_ends, chrom_ends$chr)
        fig_l <- purrr::pmap(list(sim_tb_l, pos_tb_l, chrom_ends_l),
                             .plot_regional_similarity_gg,
                     x_axis_breaks = x_axis_breaks, point_size = point_size, window_size = region_cossim@window_size,
                     stepsize = region_cossim@stepsize, oligo_correction = oligo_correction, max_cossim = max_cossim,
                     title = title, plot_rug = plot_rug)
        return(fig_l)
    }
}

#' Determine x axis breaks
#'
#' @param region_cossim A region_cossim object.
#' @param per_chrom Boolean. Determines whether to create a separate plot per chromosome.
#'
#' @return A vector with x-axis breaks.
#' @noRd
#'
.get_x_axis_breaks <- function(region_cossim, per_chrom){
    max_chr_length <- max(region_cossim@chr_lengths / 1000000)
    
    #Set x-axis breaks
    if (per_chrom == TRUE){
        x_axis_break_length <- 10
    } else{
        x_axis_break_length <- 50
    }
    
    if (max_chr_length < x_axis_break_length){
        x_axis_breaks <- x_axis_break_length
    } else{
        x_axis_breaks <- seq(x_axis_break_length, max_chr_length, by = x_axis_break_length)
    }
    return(x_axis_breaks)
}


#' Determine chromosome start and end positions.
#'
#' @param region_cossim A region_cossim object.
#'
#' @return A tibble containing the start and end position of chromosomes.
#' @noRd
#' @importFrom magrittr %>%
#'
.calc_chrom_ends <- function(region_cossim){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    chr <- NULL
    
    
    chr_pos_tb <- tibble::tibble("chr" = names(region_cossim@chr_lengths), "start" = 1, "end" = as.double(region_cossim@chr_lengths)) %>% 
        dplyr::mutate(chr = gsub("chr|chromosome_|chromosome|group|group_|chrom", "", chr), chr = factor(chr, levels = chr)) %>%
        tidyr::gather(key = "side", value = "pos", -chr) %>% 
        dplyr::mutate(pos_mb = pos / 1000000)
    return(chr_pos_tb)
}

#' Create a single regional similarity plot
#'
#' @param sim_tb A tibble containing the calculated similarities of the windows.
#' @param pos_tb A tibble containing the mutation positions.
#' @param chrom_ends A tibble containing the start and end position of chromosomes.
#' @param x_axis_breaks A vector with x-axis breaks.
#' @param point_size The point size used for plotting windows
#' @param window_size The number of mutations in a window.
#' @param stepsize The number of mutations that a window slides in each step.
#' @param oligo_correction Boolean describing whether the oligonucleotide
#'   frequency corrected cosine similarities should be plotted. If no correction
#'   has been applied then the regular cosine similarities will be plotted.
#' @param max_cossim Maximum cosine similarity for a window to be considered
#'   an outlier. Any window with a lower cosine similarity is given a different
#'   color.
#' @param title Optional plot title. (Default: NA). When the default option is
#'   used, the number of mutations per window and the step size are shown.
#' @param plot_rug Add a bottom rug to the plot, depicting the location of
#'   the mutations. (Default: FALSE)
#'
#' @return ggplot2 object
#' @noRd
#' @import ggplot2
#' @importFrom magrittr %>%
#'
.plot_regional_similarity_gg <- function(sim_tb, 
                                      pos_tb, 
                                      chrom_ends, 
                                      x_axis_breaks, 
                                      point_size, 
                                      window_size, 
                                      stepsize, 
                                      oligo_correction, 
                                      max_cossim,
                                      title,
                                      plot_rug){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    window_sizes_mb <- window_pos_mb <- pos_mb <- colour <- NULL
    
    
    # Create title
    if (is.na(title)){
        title <- paste0(" Mutations per window: ", window_size, "; Step size: ", stepsize)
    }
    
    #Set max y axis and other details
    if (oligo_correction){
        y_max_databased <- 1.1 * max(sim_tb$corrected_cossim)
        y_max <- max(y_max_databased, 1)
        y_lab <- "Corrected cosine similarity"
        y_min <- 0
        sim_type <- "corrected_cossim"
    } else{
        y_max <- 1
        y_lab <- "Cosine similarity"
        y_min <- 0
        sim_type <- "cossim"
    }
    
    sim_type <- sym(sim_type)
    
    if (!.is_na(max_cossim)){
        sim_tb <- sim_tb %>% dplyr::mutate(colour = !!sim_type <= max_cossim)
        colour_lab <- "Outlier"
        colour_scale <- scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"))
    } else{
        colour_lab <- "Window size (mb)"
        
        # Limits for colours
        sim_tb <- dplyr::mutate(sim_tb, colour = ifelse(window_sizes_mb > 4, 4, window_sizes_mb)) 
        
        # Use the blues gradient from RColorBrewer. 
        # The lightest 3 colors are removed, because they were hard to see.
        colour_scale <- scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(9, "Blues")[4:9]), 
                                              limits = c(0,4), 
                                              labels = c(0, 1, 2, 3, ">4"), 
                                              breaks = c(0, 1, 2,3, 4)) 
    }
    
    # Create main figure
    fig <- ggplot(sim_tb, aes(x = window_pos_mb, y = !!sim_type)) +
        geom_blank(data = chrom_ends, aes(x = pos_mb, y = 0)) +
        geom_point(aes(colour = colour), size = point_size) +
        facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
        labs(x = "Coordinate (mb)", y = y_lab, colour = colour_lab, title = title) +
        colour_scale +
        coord_cartesian(ylim = c(y_min, y_max), expand = FALSE) +
        theme_bw() +
        scale_x_continuous(breaks = x_axis_breaks) +
        theme(text = element_text(size = 16), axis.text.x = element_text(size = 12), legend.position = "bottom")
    
    # Optionally add rug containing mutation positions.
    if (plot_rug){
        fig <- fig +
            geom_rug(data = pos_tb, aes(x = pos_mb, y = NULL), sides = "b", size = 0.005)
    }
    
    return(fig)
}
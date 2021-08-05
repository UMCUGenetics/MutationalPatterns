#' Plot regional similarity
#' 
#' Plot the cosine similarity of the mutation profiles of small genomic windows
#' with the rest of the genome.
#' 
#' @details
#' Each dot shows the cosine similarity between the mutation profiles of a
#' single window and the rest of the genome. The dots are colored based on the
#' sizes in mega bases of the windows. This size is the distance between the
#' first and last mutations of a window. The locations of the mutations are
#' plotted on the bottom of the figure. The cosine similarity can be plotted
#' both with and without oligonucleotide frequency correction. This can be done
#' for all chromosomes at once or separate plots can be made per chromosome.
#'
#' @param region_cossim A region_cossim object.
#' @param per_chrom Boolean. Determines whether to create a separate plot per chromosome. (Default: FALSE)
#' @param tri_correction Boolean describing whether the oligonucleotide
#'   frequency corrected cosine similarities should be plotted. If no correction
#'   has been applied then the regular cosine similarities will be plotted.
#'   (Default: TRUE)
#' @param max_cossim Maximum cosine similarity for a window to be considered
#'   an outlier. Any window with a lower cosine similarity is given a different
#'   color. (Default: NA)
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
#' plot_regional_similarity(regional_sims, max_cossim = 0.5)
#' 
#' ## Plot samples per chromosome
#' fig_l = plot_regional_similarity(regional_sims, per_chrom = TRUE)
#' 
plot_regional_similarity <- function(region_cossim, 
                                    per_chrom = FALSE, 
                                    tri_correction = TRUE, 
                                    max_cossim = NA){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    chr <- sims <- metadata <- NULL
    
    #Retrieve chromosome lengths.
    max_chr_length <- max(region_cossim@chr_lengths / 1000000)
    
    #Set x-axis breaks
    x_axis_breaks <- .get_x_axis_breaks(region_cossim, per_chrom)
    
    #Get pos_tb and sim_tb.
    pos_tb <- region_cossim@pos_tb
    sim_tb <- get_sim_tb(region_cossim)
    
    # Set tri_correction to FALSE if it has not been performed
    if (!("corrected_cossim" %in% colnames(sim_tb))){
        tri_correction = FALSE
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
                                    region_cossim@stepsize, tri_correction, max_cossim = max_cossim)
        return(fig)   
    } else{
        sim_tb_l <- split(sim_tb, sim_tb$chr)
        pos_tb_l <- split(pos_tb, pos_tb$chr)
        chrom_ends_l <- split(chrom_ends, chrom_ends$chr)
        fig_l <- purrr::pmap(list(sim_tb_l, pos_tb_l, chrom_ends_l),
                             .plot_regional_similarity_gg,
                     x_axis_breaks = x_axis_breaks, point_size = point_size, window_size = region_cossim@window_size,
                     stepsize = region_cossim@stepsize, tri_correction = tri_correction, max_cossim = max_cossim)
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
#' @param tri_correction Boolean describing whether the oligonucleotide
#'   frequency corrected cosine similarities should be plotted. If no correction
#'   has been applied then the regular cosine similarities will be plotted.
#' @param max_cossim Maximum cosine similarity for a window to be considered
#'   an outlier. Any window with a lower cosine similarity is given a different
#'   color.
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
                                      tri_correction, 
                                      max_cossim){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    window_sizes_mb <- window_pos_mb <- pos_mb <- colour <- NULL
    
    
    # Create title
    title <- paste0(window_size, " muts per window; stepsize: ", stepsize)
    
    #Set max y axis and other details
    if (tri_correction){
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
        sim_tb <- sim_tb %>% dplyr::mutate(colour = ifelse(window_sizes_mb > 4, 4, window_sizes_mb)) #Limits for colours
        colour_scale <- scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(9, "Blues")[4:9]), 
                                              limits = c(0,4), 
                                              labels = c(0, 1, 2, 3, ">4"), 
                                              breaks = c(0, 1, 2,3, 4)) #Use the blues gradient from RColorBrewer. The lightest 3 colors are removed, because they were hard to see.
    }
    
    fig <- ggplot(sim_tb, aes(x = window_pos_mb, y = !!sim_type)) +
        geom_blank(data = chrom_ends, aes(x = pos_mb, y = 0)) +
        geom_point(aes(colour = colour), size = point_size) +
        geom_rug(data = pos_tb, aes(x = pos_mb, y = NULL), sides = "b", size = 0.005) +
        facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
        labs(x = "Coordinate (mb)", y = y_lab, colour = colour_lab, title = title) +
        colour_scale +
        coord_cartesian(ylim = c(y_min, y_max), expand = FALSE) +
        theme_bw() +
        scale_x_continuous(breaks = x_axis_breaks) +
        theme(text = element_text(size = 16), axis.text.x = element_text(size = 12), legend.position = "bottom")
    return(fig)
}

# #Plot an additional variable in the same way as the main cossim plot.
# #The resulting graph can be combined with the main cossim plot.
# #You can use a variable from a seperate df or a column from your original gr.
# plot_sims_additional = function(sims, column = NA, df = NA, per_chrom = F){
#     
#     #Check arguments
#     check_sims(sims)
#     check_logical(per_chrom)
#     
#     #Preprocessing
#     x_axis_breaks = get_x_axis_breaks(sims, per_chrom)
#     
#     if (!is_na(df)){
#         variable_tb = transform_to_variable_tb(sims, df, column)
#     } else if (!is_na(column)){
#         variable_tb = get_variable_tb(sims, column)
#     } else{
#         stop("Please provide either a column in the original gr that can be used\n
#              or a seperate df.", call. = F)
#     }
#     
#     chrom_ends = calc_chrom_ends(sims)
#     
#     #Plot either discrete or continous depending on user input.
#     if ("variable" %in% colnames(variable_tb)){
#         plot_sims_additional_continuous(variable_tb, per_chrom, x_axis_breaks, column, chrom_ends)
#     } else{
#         plot_sims_additional_discrete(variable_tb, per_chrom, x_axis_breaks, chrom_ends)
#     }
# }
# 
# 
# 
# #Helper function for plot_sims_additional. Converts a df to the right format.
# transform_to_variable_tb = function(sims, df, column){
#     
#     if (!"chr" %in% colnames(df)){
#         stop("Please provide a data frame with a column named `chr`.", call. = F)
#     }
#     
#     variable_tb = df %>% 
#         as_tibble() %>%
#         dplyr::filter(chr %in% levels(sims@pos_tb$chr))
#     if (nrow(variable_tb) == 0){
#         stop("The supplied df has no seqlevels in common with the sims object.\n
#              Please ensure they use the same seqlevelstyle.", .call = F)
#     }
#     
#     if ("start" %in% colnames(variable_tb)){#Use start and end
#         variable_tb = variable_tb %>% 
#             dplyr::mutate(chr = factor(chr, levels = levels(sims@pos_tb$chr)),
#                           start_mb = start / 1000000,
#                           end_mb = end / 1000000)
#         levels(variable_tb$chr) = gsub("chr|Chr|CHR", "", levels(variable_tb$chr))
#         
#         if (is.na(column)){
#             variable_tb = variable_tb %>% 
#                 dplyr::select(chr, start, end, start_mb, end_mb)
#         } else{
#             variable_tb = variable_tb %>% 
#                 dplyr::select(chr, start, end, start_mb, end_mb, variable = column)
#         }
#         
#         return(variable_tb)
#     } else{#Use position
#         variable_tb = variable_tb %>% 
#             dplyr::mutate(chr = factor(chr, levels = levels(sims@pos_tb$chr)),
#                           pos_mb = pos / 1000000)
#         levels(variable_tb$chr) = gsub("chr|Chr|CHR", "", levels(variable_tb$chr))
#         
#         if (is.na(column)){
#             variable_tb = variable_tb %>% 
#                 dplyr::select(chr, pos, pos_mb)
#         } else{
#             variable_tb = variable_tb %>% 
#                 dplyr::select(chr, pos, pos_mb, variable = column)
#         }
#         
#         return(variable_tb)
#     }
# }
# 
# #Helper function for plot_sims_additional. Get the correct collumn from the gr and converts it to the right format.
# get_variable_tb = function(sims, column){
#     gr = sims@gr
#     variable_tb = tibble("chr" = as.vector(seqnames(gr)), "pos" = start(gr), "variable" = mcols(gr)[,column]) %>%
#         dplyr::mutate(chr = factor(chr, levels = levels(sims@pos_tb$chr)),
#                       pos_mb = pos / 1000000)
#     levels(variable_tb$chr) = gsub("chr|Chr|CHR", "", levels(variable_tb$chr))
#     return(variable_tb)
# }
# 
# 
# #Helper function for plot_sims_additional. Works on continuous values.
# plot_sims_additional_continuous = function(variable_tb, per_chrom, x_axis_breaks, column, chrom_ends){
#     
#     if (per_chrom == F){
#         fig = plot_sims_additional_continuous_gg(variable_tb, x_axis_breaks, column, chrom_ends)
#         return(fig)
#     } else if (per_chrom == T){
#         variable_tb_l = split(variable_tb, variable_tb$chr)
#         chrom_ends_l = split(chrom_ends, chrom_ends$chr)
#         fig_l = purrr::map2(variable_tb_l, chrom_ends_l, function(x, y) plot_sims_additional_continuous_gg(x, x_axis_breaks, column, y))
#         return(fig_l)
#     }
# }
# 
# #Helper function for plot_sims_additional_continuous. Performs the actual plotting
# plot_sims_additional_continuous_gg = function(variable_tb, x_axis_breaks, column, chrom_ends){
#     
#     #Check if data is empty
#     if (nrow(variable_tb) > 0){
#         data_geom = geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))
# 
#     } else{
#         data_geom = geom_blank(data = chrom_ends)
#     }
#     
#     
#     fig = ggplot(variable_tb, aes(x = pos_mb, y = variable)) +
#         geom_blank(data = chrom_ends, aes(x = pos_mb, y = 0)) +
#         data_geom +
#         facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
#         theme_bw() +
#         labs(x = "Coordinate (mb)", y = column) +
#         scale_x_continuous(breaks = x_axis_breaks) +
#         theme(text = element_text(size = 18), axis.text.x = element_text(size = 12), legend.position = "bottom")
#     return(fig)
# }
# 
# #Helper function for plot_sims_additional. Works on discrete values.
# plot_sims_additional_discrete = function(variable_tb, per_chrom, x_axis_breaks, chrom_ends){
#     if (per_chrom == F){
#         fig = plot_sims_additional_discrete_gg(variable_tb, x_axis_breaks, chrom_ends)
#         return(fig)
#     } else if (per_chrom == T){
#         variable_tb_l = split(variable_tb, variable_tb$chr)
#         chrom_ends_l = split(chrom_ends, chrom_ends$chr)
#         fig_l = purrr::map2(variable_tb_l, chrom_ends_l, function(x, y) plot_sims_additional_discrete_gg(x, x_axis_breaks, y))
#         return(fig_l)
#     }
# }
# 
# 
# #Helper function for plot_sims_additional_discrete. Performs the actual plotting
# plot_sims_additional_discrete_gg = function(variable_tb, x_axis_breaks, chrom_ends){
#     
#     #Check if data is empty and set a data geom
#     if (nrow(variable_tb) == 0){
#         data_geom = geom_blank(data = chrom_ends)
#     } else if ("start" %in% colnames(variable_tb)){
#         data_geom = geom_rect(aes(xmin = start_mb, xmax = end_mb))
#     } else{
#         data_geom = geom_linerange(aes(x = pos_mb), size = 0.5)
#         
#     }
#     
#     
#     fig = ggplot(variable_tb, aes(ymin = 0, ymax = 1)) +
#         geom_blank(data = chrom_ends, aes(x = pos_mb)) +
#         data_geom +
#         facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
#         theme_bw() +
#         labs(x = "Coordinate (mb)", y = "") +
#         scale_x_continuous(breaks = x_axis_breaks) +
#         theme(text = element_text(size = 18), axis.text.x = element_text(size = 12), legend.position = "bottom", 
#               axis.text.y = element_blank(), axis.ticks.y = element_blank())
#     return(fig)
#     
# }
# 
# #Function to combine the main cossim figure with a figure showing an additional variable.
# arrange_sims_other_figs = function(sims_fig, other_fig){
#     
#     check_ggplot(sims_fig)
#     check_ggplot(other_fig)
# 
#     arrangeable_sims_fig = sims_fig +
#         theme(axis.title.x = element_blank(),
#               axis.text.x = element_blank(),
#               axis.ticks.x = element_blank()
#         )
#     arrangeable_other_fig = other_fig +
#         theme(strip.background = element_blank(),
#               strip.text.x = element_blank())
#     fig = ggarrange(arrangeable_sims_fig,
#                     arrangeable_other_fig,
#                     nrow = 2,
#                     align = "v",
#                     heights = c(2, 0.7),
#                     legend = "right",
#                     common.legend = T
#     )
#     return(fig)
# }
# 
# 
# #Plot sims without sample. Works on the st. jude data.
# write_plot_sims_without_sample = function(sims, sample){
#     gr = sims@gr
#     gr_nosample = gr[gr$sampleID != sample]
#     sims_nosample = get_window_cossim(gr_nosample, window_size = sims@window_size, stepsize = sims@stepsize, max_window_size_gen = sims@max_window_size_gen, ref_genome = sims@ref_genome, chromosomes = sims@chromosomes, tri_correction = F)
#     sims_nosample_fig = plot(sims_nosample)
#     sims_nosample_fig_l = plot(sims_nosample, per_chrom = T)
#     
#     pdf(paste0("consistent_lowsim_regions/local_mutpattern_no", sample, ".pdf"), width = 25, useDingbats = F)
#     print(sims_nosample_fig)
#     print(sims_nosample_fig_l)
#     dev.off()
#     invisible(0)
# }

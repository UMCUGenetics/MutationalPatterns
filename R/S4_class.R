#S4 class specifications and methods


#' @importFrom methods setClass setGeneric setMethod
NULL

#' An S4 class to store the results of a regional mutation pattern similarity
#' analysis
#'
#' @slot sim_tb A tibble containing the calculated similarities of the windows.
#' @slot pos_tb A tibble containing the mutation positions.
#' @slot chr_lengths Vector containing the chromosome lengths.
#' @slot window_size The number of mutations in a window.
#' @slot max_window_size_gen The maximum size of a window before it is removed.
#' @slot ref_genome BSgenome reference genome object
#' @slot muts_per_chr Vector containing the number of mutations per chromosome.
#' @slot mean_window_size The mean length of the genome covered by the windows.
#' @slot stepsize The number of mutations that a window slides in each step.
#' @slot extension The number of bases, that's extracted upstream and
#'   downstream of the base substitutions, to create the mutation matrices.
#' @slot chromosomes Vector of chromosome/contig names of the reference genome
#'   to be plotted.
#' @slot exclude_self_mut_mat Boolean describing whether the mutations in a
#'   window should be subtracted from the global mutation matrix.
#'
#'
#' @export
setClass("region_cossim", slots = c("sim_tb", "pos_tb", "chr_lengths", "window_size", 
                                   "max_window_size_gen", "ref_genome", "muts_per_chr", 
                                   "mean_window_size", "stepsize", "extension", "chromosomes", "exclude_self_mut_mat"))

#' An S4 method to show an instance of the region_cossim class.
#' 
#' @param object A region_cossim object.
#' 
#' @export
setMethod("show", signature = "region_cossim", function(object){
    cat(is(object), " object defined on genome: ", S4Vectors::metadata(object@ref_genome)$genome, "\n",
        "Mutations per window: ", object@window_size, "\n",
        "Step size: ", object@stepsize, "\n",
        "Max allowed window size in bp: ", object@max_window_size_gen, "\n",
        "Mean window size: ", object@mean_window_size, "\n",
        "Mutations per chromosome: \n", sep = "")
    cat(object@muts_per_chr)
})

#' An S4 generic to get the sim_tb from a region_cossim object.
#' 
#' @param x A region_cossim object
#'
#' @return A tibble containing the calculated similarities of the windows.
#' 
#' @export
setGeneric("get_sim_tb", function(x) standardGeneric("get_sim_tb"))


#' An S4 method for the get_sim_tb generic
#' 
#' @param region_cossim A region_cossim object
#' 
#' @return A tibble containing the calculated similarities of the windows.
#'
#' @export
#' @describeIn get_sim_tb Get the sim_tb from a region_cossim object.
setMethod("get_sim_tb", "region_cossim", function(x) x@sim_tb)



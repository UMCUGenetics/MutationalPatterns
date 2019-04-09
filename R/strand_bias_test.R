#' Significance test for strand asymmetry
#'
#' This function performs a Poisson test for the ratio between mutations on 
#' each strand
#' 
#' @param strand_occurrences Dataframe with mutation count per strand, result
#' from strand_occurrences()
#' @return Dataframe with poisson test P value for the ratio between the
#' two strands per group per base substitution type.
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#' @importFrom plyr .
#' @importFrom plyr ddply
#' @importFrom plyr summarise
#' @importFrom stats "poisson.test"
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## following mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' tissue <- c("colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver")
#'
#' ## Perform the strand bias test.
#' strand_counts = strand_occurrences(mut_mat_s, by=tissue)
#' strand_bias = strand_bias_test(strand_counts)
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurrences}},
#' \code{\link{plot_strand_bias}}
#'
#' @export

strand_bias_test = function(strand_occurrences)
{
    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    group = NULL
    type = NULL
    strand = NULL
    variable = NULL

    # statistical test for strand ratio
    # poisson test
    df_strand = reshape2::dcast(melt(strand_occurrences),
                                group + type ~ strand,
                                sum,
                                subset = plyr::.(variable == "no_mutations"))
    
    df_strand$total = df_strand[,3] + df_strand[,4]
    df_strand$ratio = df_strand[,3] / df_strand[,4]
    df_strand$p_poisson = apply(df_strand, 1, function(x) poisson.test(c(as.numeric(x[3]), as.numeric(x[4])), r=1)$p.value)
    df_strand$significant[df_strand$p_poisson < 0.05] = "*"
    df_strand$significant[df_strand$p_poisson >= 0.05] = " "

    return(df_strand)
}

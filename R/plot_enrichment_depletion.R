#' Plot enrichment/depletion of mutations in genomic regions
#' 
#' @param df Dataframe result from enrichment_depletion_test()
#' @return Plot with two parts. 1: Barplot with no. mutations expected and
#' observed per region. 2: Effect size of enrichment/depletion
#' (log2ratio) with results significance test.
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom reshape2 melt
#'
#' @examples
#' ## See the 'genomic_distribution()' example for how we obtained the
#' ## following data:
#' distr <- readRDS(system.file("states/distr_data.rds",
#'                     package="MutationalPatterns"))
#' 
#' tissue = c( "colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver" )
#'
#' ## Perform the enrichment/depletion test.
#' distr_test = enrichment_depletion_test(distr, by = tissue)
#' distr_test2 = enrichment_depletion_test(distr)
#'
#' ## Plot the enrichment/depletion
#' plot_enrichment_depletion(distr_test)
#' plot_enrichment_depletion(distr_test2)
#'
#' @seealso
#' \code{\link{enrichment_depletion_test}},
#' \code{\link{genomic_distribution}}
#'
#' @export

plot_enrichment_depletion = function(df)
{
    df2 = melt(df[,c(1,2,6,8)], id = c("by", "region"))

    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    value = NULL
    variable = NULL
    observed = NULL
    expected = NULL
    significant = NULL

    # Part 1: No. mutations expected and observed per region
    plot1 =  ggplot(df2, aes(x=by,
                                y=value,
                                fill=by,
                                group=variable,
                                alpha=variable)) +
        geom_bar(colour="black",
                    stat="identity",
                    position=position_dodge()) +
        facet_grid(~ region) +
        theme_bw()  +
        theme(axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                legend.title=element_blank()) +
        xlab("") +
        ylab("No. mutations") +
    scale_x_discrete(breaks=NULL)

    # determine max y value for plotting
    # = log2 ratio with pseudo counts
    max = round(max(abs(log2((df$observed+0.1) / (df$expected+0.1)))),
                digits = 1) + 0.1

    # Part 2: effect size of enrichment/depletion with significance test
    plot2 = ggplot(data=df, aes(x=by,
                                y=log2((observed+0.1)/(expected+0.1)),
                                fill=by)) +
        geom_bar(colour="black",
                    stat="identity",
                    position=position_dodge()) +
        scale_y_continuous(limits=c(-max, max)) +
        geom_text(
            aes(x = by,
                y = log2((observed+0.1) / (expected+0.1)),
                label = significant,
                vjust = ifelse(sign(log2((observed+0.1) /
                                            (expected+0.1))) > 0, 0.5, 1)),
                size = 8, position = position_dodge(width = 1)) +
        facet_grid(~ region) +
        theme_bw() +
        theme(axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                legend.title = element_blank()) +
        xlab("") +
        ylab("log2(observed/expected)") +
        scale_x_discrete(breaks = NULL)

    output <- cowplot::plot_grid (plot1, plot2, ncol=1, nrow=2, rel_heights = c(2,1.2))
    return(output)
}

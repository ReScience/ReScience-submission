#' Run standard target analysis
#' 
#' Load the data from the text files. Extract and resize regions. Call the
#' target package function direct_targets 
#' 
#' @param peak_file A string. The path to the peaks file.
#' @param peak_filter A logical vector to filter peak_file. Default TRUE 
#' @param deg_file A string. The path to the differential expression file.
#' @param deg_filter A logical vector to filter deg_file. Default TRUE 
#' @param genome_file A string. The path to the genome file
#' @param genome_filter A logical vector to filter genome_file. Default TRUE 
#' @param base An integer. The length or the genomic region of interest in bp.
#' @param regions_col A string. The name of the names column in deg_file.
#' @param stats_col A string. The name of the statistics column in deg_file.
#'
#' @return A GRanges object. The objects contains the resized genomic regions
#' of interest (seqnames, ranges and strand). These are extracted in two steps.
#' First, the coordinates of all regions in deg_file are extracted from the 
#' genome file and resized to the base length flanks on either side of their
#' start. Second, only regions that intersect with the peaks in peak_file are
#' kept. The metadata part contains all the corresponding columns from deg_file
#' and the five extra columns. The five extra columns are
#' \describe{
#' \item score The transformed distance between the peak and the region start
#' \item score_rank The rank of the score
#' \item stat The statistics column from deg_file indicated by stats_col
#' \item stat_rank The rank of the stat
#' \item rank_product The product of score_rank and stat_rank
#' }
#' 
#' @examples 
#' library(tidyverse)
#' run_target(
#'     peak_file = 'data/beta_input/3656_peaks.bed',
#'     deg_file = 'data/beta_input/AR_diff_expr.xls',
#'     genome_file = 'data/beta_input/hg19.refseq',
#'     regions_col = 'refseq',
#'     stats_col = 't',
#'     base = 100000
#'     )
#'     
#' @export
run_target <- function(peak_file, peak_filter = TRUE,
                       deg_file, deg_filter = TRUE, 
                       genome_file, genome_filter = TRUE,
                       base, regions_col, stats_col,
                       n = NULL) {
    # read peaks file, and
    # transform to genomic ranges
    peaks <- readr::read_tsv(peak_file, 
                             col_names = c('seqnames', 'start', 'end', 'name', 'pval')) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # read genome file, and
    # transform to granges, and
    # extract regions, and
    # back to tbl
    genome <- readr::read_tsv(genome_file,
                              col_names = c('refseq', 'chrom', 'strand', 'start', 'end', 'symbol'),
                              skip = 1) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
        GenomicRanges::promoters(upstream = base, downstream = base) %>%
        dplyr::as_tibble()
    
    # read differential expression file, and
    # clean refseq, and
    # merge with genom, and
    # transform to granges
    regions <- readr::read_tsv(deg_file,
                               col_names = c("refseq", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val","B" ),
                               skip = 1) %>%
        dplyr::mutate(refseq = str_split(refseq, '\\_at', simplify = TRUE)[, 1])
    
    if (is.numeric(deg_filter)) {
        regions <- dplyr::bind_rows(
            dplyr::top_n(regions, deg_filter, t),
            dplyr::top_n(regions, -deg_filter, t)
        )
    }
    regions <- regions %>%
        dplyr::inner_join(genome) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # call direct tragets
    target::direct_targets(peaks = peaks,
                           regions = regions, 
                           regions_col = regions_col,
                           stats_col = stats_col,
                           base = base,
                           n = n)
}

#' Plot the ECDF of the direct targets in groups
#' 
#' Plot the ECDF of the regulatory potential of the predicted direct targets in
#' groups defined by the quantiles.
#'
#' @param res_file A string. The path to the rds file produced by run_target.
#' @param cut_quantiles A vector of numerics. The four quantiles to use in cut.
#' @param cut_top An integer to select the top n genes in both directions.
#'
#' @return A ggplot object. The graph show the ECDF of the regulatory potential
#' of three groups of regions based the logFC.
#' 
#' @examples 
#' library(tidyverse)
#' plot_target(res_file = 'data/target_output/basic_AR_3656.rds',
#'             cut_quantiles = c(0, .1, .9, 1))
#'             
#' @export
plot_target <- function(res_file, cut_quantiles = NULL, cut_top = NULL) {
    # read results file
    gr <- readr::read_rds(res_file)
    df <- dplyr::as_tibble(gr)
    
    # transform into tbl, and
    # cut the logFC column by cut quantiles, and
    # make plot
    if (is.null(cut_top)) {
        df <- df %>%
            dplyr::mutate(group = cut(logFC,
                                      breaks = stats::quantile(logFC, cut_quantiles),
                                      labels = c('Down', 'None', 'Up')))
        
    } else if (is.null(cut_quantiles)) {
        top <- dplyr::top_n(df, cut_top, logFC) %>% dplyr::mutate(group = 'Up')
        bot <- dplyr::top_n(df, -cut_top, logFC) %>% dplyr::mutate(group = 'Down')
        df <- dplyr::full_join(df, top) %>%
            dplyr::full_join(bot) %>%
            dplyr::mutate(group = ifelse(is.na(group), 'None', group))
    }
    df %>%
        ggplot2::ggplot(aes(x = score_rank, color = group)) +
        ggplot2::stat_ecdf() +
        ggplot2::theme_bw() +
        ggplot2::geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
        ggplot2::theme(panel.grid = element_blank(),
                       legend.position = c(.7,.3),
                       legend.background = element_rect(fill = NA),
                       panel.border = element_rect(size = 1.5)) +
        ggplot2::labs(x = 'Regulatory Potential', y = 'ECDF', color = '') +
        ggplot2::scale_color_manual(values = c('darkblue', 'gray', 'darkred'))
}

#' Extract top direct targets
#' 
#' Extart the top n direct target as predicted by the target package function
#' direct_targets
#' 
#' @param res_file A string. The path to the rds file produced by run_target.
#' @param n An integet.
#'
#' @return A data.frame with six columns:
#' \describe{
#' \item seqnames The chromosome name
#' \item start The start coordinate of the region
#' \item end The end coordinate of the region
#' \item refseq The REFSEQ identifier of the rgion
#' \item symbol The GENE SYMBOL of the region
#' \item rank The calculated rank of the region
#' }
#' 
#' @examples 
#' top_target(res_file = 'data/target_output/basic_AR_3656.rds',
#'            n = 3)
#'            
#' @export
top_target <- function(res_file, n) {
    # load the rds file,
    # transform to tibble,
    # arrance by rank,
    # select relevant columns, and
    # return top n
    readr::read_rds(res_file) %>%
        dplyr::as_tibble() %>%
        dplyr::arrange(rank) %>%
        dplyr::select(seqnames, start, end, refseq, symbol, rank) %>%
        head(n)
}

#' Calculate the concordance of predicted sets of targets
#' 
#' Calculate the fraction of intersections between multiple sets of predicted 
#' direct targets. It defines one set as default and compares it to all other
#' sets.
#'
#' @param res_file A string. The path to the rds file produced by run_target.
#' @param compare_to A string. The name of the set to compare the others to.
#' @param sets A vector of integers. The lengths of the sets to compare at.
#' 
#' @return A data.frame of three columns:
#' \describe{
#' \item type The name of the set
#' \item fraction The fraction of intersection in this set to the default
#' \item set The size of the set
#' }
#' 
#' @examples
#' sets <- seq(1, 10000, 99)
#' 
#' concordance('data/target_output/variants_AR_3656.rds',
#'             compare_to = 'defaults',
#'             sets = sets)
#'                  
#' @export
concordance <- function(res_file, compare_to, sets) {
    # load the rds file,
    # arrange by ranks, and
    # split in a list by variant
    dt <- readr::read_rds(res_file) %>%
        dplyr::arrange(rank) %>%
        with(split(.$symbol, .$variant))
    
    # for each variant,
    # calculate the fraction of interesection with the rest
    purrr::map_df(dt, function(x) {
        # extract the default variant
        vec1 <- dt[[compare_to]]
        
        # extract the other set
        vec2 <- x
        
        # make empty vectors
        res <- vector('numeric', length = length(sets))
        set <- vector('integer', length = length(sets))
        
        # at each set length 
        # get the fraction of the intersection
        for(i in seq_along(sets)) {
            n <- sets[i]
            
            x <- vec1[1:n]
            y <- vec2[1:n]
            
            res[i] <- length(intersect(x, y))/length(unique(x))
            set[i] <- n
        }
        
        # tidy output
        tibble(fraction = res,
               set = set)
    },
    .id = 'type')
}

#' Plot the concordance of predicted sets of targets
#' 
#' Plot the fraction of intersections between multiple sets of predicted 
#' direct targets. It defines one set as default and compares it to all other
#' sets. The calculations are made by the concordence function
#'
#' @param res_file A string. The path to the rds file produced by run_target.
#' @param compare_to A string. The name of the set to compare the others to.
#' @param sets A vector of integers. The lengths of the sets to compare at.#'
#'
#' @return A ggplot graph. The graph shows the fraction of intersection of 
#' predicted target genes between a default set given by compare to and the 
#' rest at set sizes given by sets.
#' 
#' @examples
#' sets <- seq(1, 10000, 99)
#' 
#' plot_concordance('data/target_output/variants_AR_3656.rds',
#'                  compare_to = 'defaults',
#'                  sets = sets)
#'   
#' @export
plot_concordance <- function(res_file, compare_to, sets) {
    # calculate the concordence
    df <- concordance(res_file,
                      compare_to,
                      sets)
    
    # make a graph
    ggplot2::ggplot(df, aes(x = set, y = fraction, color = type)) +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = element_blank(),
                       legend.background = element_rect(fill = NA),
                       legend.position = c(.75,.3),
                       legend.text = element_text(size = 6),
                       panel.border = element_rect(size = 1.5)) +
        ggplot2::labs(x = 'Number of Genes', y = 'Fraction', color = '')
}

#' Test predicted factor functions at diffenent cut-offs
#' 
#' Test the significance of the predicted functions of the factors at different
#' cut-offs for the grouping of targets.
#'
#' @param res_file A string. The path to the rds file produced by run_target.
#' @param quants A list. Each item of the list is a vector of numerics. The 
#' four quantiles to use in cut.
#'
#' @return A data.frame with six columns:
#' \describe{
#' \item quantiles The names of the quants list items
#' \item direction The direction of the regulation; down, none or up
#' \item statistic, p.value, method and alternative As in test_predicitons
#' }
#'
#' @examples
#' quants <- list(
#'     ten = c(0, .1, .9, 1),
#'     twenty = c(0, .2, .8, 1),
#'     thirty = c(0, .3, .7, 1),
#'     fourty = c(0, .4, .6, 1)
#' )
#' 
#' test_cutoffs('data/target_output/basic_AR_3656.rds',
#'              quants = quants)
#' 
#' @export
test_cutoffs <- function(res_file, quants) {
    # load the rds file,
    # transform to tibble, and
    # remove missing values
    dt <- readr::read_rds(res_file) %>% 
        dplyr::as_tibble() %>%
        stats::na.omit()
    
    # for each set of quantiles,
    # group targets by quantiles,
    # run test predictions, and
    # tidy
    purrr::map_df(quants, function(x) {
        group <- cut(dt$stat,
                     breaks = stats::quantile(dt$stat, x),
                     labels = c('down', 'none', 'up'))
        
        map_df(c(up = 'up', down = 'down'), function(y) {
            target::test_predictions(dt$score_rank,
                                     group = group,
                                     compare = c(y, 'none')) %>%
                broom::tidy()
        }, .id = 'direction')
    }, .id = 'quantiles')
}
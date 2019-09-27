run_target <- function(peak_file, deg_file, genome_file, base, regions_col, stats_col) {
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
        dplyr::mutate(refseq = str_split(refseq, '\\_at', simplify = TRUE)[, 1]) %>%
        dplyr::inner_join(genome) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # call direct tragets
    target::direct_targets(peaks = peaks,
                           regions = regions, 
                           regions_col = regions_col,
                           stats_col = stats_col,
                           base = base)
}

plot_target <- function(res_file, cut_quantiles) {
    # read results file
    gr <- read_rds(res_file)
    
    # transform into tbl, and
    # cut the logFC column by cut quantiles, and
    # make plot
    dplyr::as_tibble(gr) %>%
        dplyr::mutate(group = cut(logFC,
                           breaks = stats::quantile(logFC, cut_quantiles),
                           labels = c('Down', 'None', 'Up'))) %>%
        stats::na.omit() %>%
        ggplot2::ggplot(aes(x = score_rank, color = group)) +
        ggplot2::stat_ecdf() +
        ggplot2::theme_bw() +
        ggplot2::geom_hline(yintercept = c(0, 1), lty = 2, color = 'gray') +
        ggplot2::theme(panel.grid = element_blank(),
                       legend.position = c(.7,.3),
                       legend.background = element_rect(fill = NA),
                       panel.border = element_rect(size = 1.5)) +
        ggplot2::labs(x = 'Regulatory Potential', y = 'ECDF', color = '') +
        ggplot2::scale_color_manual(values = c('darkgreen', 'gray', 'darkred'))
}

top_target <- function(res_file, n) {
    readr::read_rds(res_file) %>%
        dplyr::as_tibble() %>%
        dplyr::arrange(rank) %>%
        dplyr::select(seqnames, start, end, refseq, symbol, rank) %>%
        head(n)
}

concordance <- function(res_file, compare_to, sets) {
    dt <- readr::read_rds(res_file) %>%
        dplyr::arrange(rank) %>%
        with(split(.$symbol, .$variant))
    
    purrr::map_df(dt, function(x) {
        vec1 <- dt[[compare_to]]
        vec2 <- x
        res <- vector('numeric', length = length(sets))
        set <- vector('integer', length = length(sets))
        for(i in seq_along(sets)) {
            n <- sets[i]
            
            x <- vec1[1:n]
            y <- vec2[1:n]
            
            res[i] <- length(intersect(x, y))/length(unique(x))
            set[i] <- n
        }
        tibble(fraction = res,
               set = set)
    },
    .id = 'type')
}

plot_concordance <- function(res_file, compare_to, sets) {
    df <- concordance(res_file, compare_to, sets)
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

test_cutoffs <- function(res_file, quants) {
    dt <- readr::read_rds(res_file) %>% 
        dplyr::as_tibble() %>%
        stats::na.omit()
    
    
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
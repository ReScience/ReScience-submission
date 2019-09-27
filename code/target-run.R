# load required libraries ----
library(tidyverse)
library(GenomicRanges)
library(target)

source('code/function-def.R')

# make directories
if(!dir.exists('data/')) dir.create('data/')
if(!dir.exists('data/target_output/')) dir.create('data/target_output/')

# call target
run_target(
    peak_file = 'data/beta_input/3656_peaks.bed',
    deg_file = 'data/beta_input/AR_diff_expr.xls',
    genome_file = 'data/beta_input/hg19.refseq',
    regions_col = 'refseq',
    stats_col = 't',
    base = 100000
) %>%
    write_rds('data/target_output/basic_AR_3656.rds')

run_target(
    peak_file = 'data/beta_input/349_peaks.bed',
    deg_file = 'data/beta_input/ESR1_diff_expr.xls',
    genome_file = 'data/beta_input/hg19.refseq',
    regions_col = 'refseq',
    stats_col = 't',
    base = 100000
) %>%
    write_rds('data/target_output/basic_ESR1_349.rds')

run_target(
    peak_file = 'data/beta_input/5795_peaks.bed',
    deg_file = 'data/beta_input/Tet1_diff_expr.xls',
    genome_file = 'data/beta_input/mm9.refseq',
    regions_col = 'refseq',
    stats_col = 't',
    base = 100000
) %>%
    write_rds('data/target_output/basic_Tet1_5795.rds')

# call target with varying parameters
param <- list(
    `defaults` = list(
        genome_file = 'data/beta_input/hg19.refseq',
        regions_col = 'refseq',
        stats_col = 't',
        base = 100000
    ),
    hg18 = list(
        genome_file = 'data/beta_input/hg18.refseq',
        regions_col = 'refseq',
        stats_col = 't',
        base = 100000
    ),
    `fold-change` = list(
        genome_file = 'data/beta_input/hg19.refseq',
        regions_col = 'refseq',
        stats_col = 'logFC',
        base = 100000
    ),
    `50kb` = list(
        genome_file = 'data/beta_input/hg19.refseq',
        regions_col = 'refseq',
        stats_col = 't',
        base = 50000
    )
)

map_dfr(param, 
    function(x) {
        run_target(
            peak_file = 'data/beta_input/3656_peaks.bed',
            deg_file = 'data/beta_input/AR_diff_expr.xls',
            genome_file = x$genome_file,
            regions_col = x$regions_col,
            stats_col = x$stats_col,
            base = x$base
        ) %>%
            as_tibble()
    },
    .id = 'variant') %>%
    write_rds('data/target_output/variants_AR_3656.rds')

map_dfr(param, 
        function(x) {
            run_target(
                peak_file = 'data/beta_input/349_peaks.bed',
                deg_file = 'data/beta_input/ESR1_diff_expr.xls',
                genome_file = x$genome_file,
                regions_col = x$regions_col,
                stats_col = x$stats_col,
                base = x$base
            ) %>%
                as_tibble()
        },
        .id = 'variant') %>%
    write_rds('data/target_output/variants_ESR1_349.rds')

param_mouse <- list(
    defaults = list(
        genome_file = 'data/beta_input/mm9.refseq',
        regions_col = 'refseq',
        stats_col = 't',
        base = 100000
    ),
    mm10 = list(
        genome_file = 'data/beta_input/mm10.refseq',
        regions_col = 'refseq',
        stats_col = 't',
        base = 100000
    ),
    `fold-change` = list(
        genome_file = 'data/beta_input/mm9.refseq',
        regions_col = 'refseq',
        stats_col = 'logFC',
        base = 100000
    ),
    `50kb` = list(
        genome_file = 'data/beta_input/mm9.refseq',
        regions_col = 'refseq',
        stats_col = 't',
        base = 50000
    )
)

map_dfr(param_mouse, 
        function(x) {
            run_target(
                peak_file = 'data/beta_input/5795_peaks.bed',
                deg_file = 'data/beta_input/Tet1_diff_expr.xls',
                genome_file = x$genome_file,
                regions_col = x$regions_col,
                stats_col = x$stats_col,
                base = x$base
            ) %>%
                as_tibble()
        },
        .id = 'variant') %>%
    write_rds('data/target_output/variants_Tet1_5795.rds')

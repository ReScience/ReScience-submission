# load required libraries
library(tidyverse)
library(xtable)

source('code/function-def.R')

# make directories
if(!dir.exists('article/')) dir.create('article/')
if(!dir.exists('article/tables/')) dir.create('article/tables/')

# make table

list(AR = 'data/target_output/basic_AR_3656.rds',
     ESR1 = 'data/target_output/basic_ESR1_349.rds',
     Tet1 = 'data/target_output/basic_Tet1_5795.rds') %>%
    map_dfr(function(x) {
        top_target(x, n = 3)
    }, 
    .id = 'factor') %>%
    mutate(refseq = str_replace(refseq, '\\_', '\\\\_'),
           factor = ifelse(duplicated(factor), '', factor),
           rank = format(rank, scientific = TRUE, digits = 3)) %>%
    setNames(c('Factor', 'Chr', 'Start', 'End', 'Refseq', 'Symbol', 'Rank')) %>%
    xtable(align = 'clllllll') %>%
    print(floating = FALSE,
          include.rownames = FALSE,
          booktabs = TRUE,
          sanitize.text.function = identity,
          comment = FALSE,
          math.style.exponents = TRUE,
          add.to.row = list(pos = list(3, 6), command = rep('\\midrule ', 2)),
          file = 'article/tables/heads.tex')

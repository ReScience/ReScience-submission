# load libraries
library(tidyverse)
library(xtable)

source('code/function-def.R')

quants <- list(
    ten = c(0, .1, .9, 1),
    twenty = c(0, .2, .8, 1),
    thirty = c(0, .3, .7, 1),
    fourty = c(0, .4, .6, 1)
)

header <- paste0(
    "\\multirow{2}{*}{Factor} & \\multirow{2}{*}{Quantiles} &",
    "\\multicolumn{2}{c}{Down} & \\multicolumn{2}{c}{Up}\\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
    "&& stat & p-value & stat & p-value \\\\"
)

list(AR = 'data/target_output/basic_AR_3656.rds',
     ESR1 = 'data/target_output/basic_ESR1_349.rds',
     Tet1 = 'data/target_output/basic_Tet1_5795.rds') %>%
    map_df(~test_cutoffs(.x, quants),
           .id = 'factor') %>%
    select(-method, -alternative) %>%
    mutate(statistic = round(statistic, 2),
           p.value = format(p.value, scientific = TRUE, digits = 1)) %>%
    unite('value', statistic, p.value) %>%
    spread(direction, value) %>%
    separate(down, c('down_stat', 'down_pval'), sep = '_') %>%
    separate(up, c('up_stat', 'up_pval'), sep = '_') %>%
    mutate(quantiles = factor(quantiles, levels = c('ten', 'twenty', 'thirty', 'fourty'))) %>%
    arrange(factor, quantiles) %>%
    mutate(factor = ifelse(duplicated(factor), '', factor),
           quantiles = case_when(quantiles == 'ten' ~ '(0.1-0.9)',
                                 quantiles == 'twenty' ~ '(0.2-0.8)',
                                 quantiles == 'thirty' ~ '(0.3-0.7)',
                                 quantiles == 'fourty' ~ '(0.4-0.6)')) %>%
    xtable(align = 'cllcccc') %>%
    print(floating = FALSE,
          include.rownames = FALSE,
          include.colnames = FALSE,
          booktabs = TRUE,
          sanitize.text.function = identity,
          comment = FALSE,
          math.style.exponents = TRUE,
          add.to.row = list(pos = list(0, 4, 8), command = c(header, rep('\\midrule ', 2))),
          file = 'article/tables/cutoffs.tex')
    

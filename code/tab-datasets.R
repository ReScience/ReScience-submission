# load required libraries
library(tidyverse)
library(xtable)

# make table

tibble(
    Factor = c('AR', 'ESR1', 'Tet1'),
    `Cell Line` = c('LNCaP', 'MCF-7', 'ES Cells'),
    Genome = c('hg19', 'hg19', 'mm9'),
    Treatment = c('DHT', 'E2', 'Knockdown'),
    `Binding Data` = c('\\cite{Wang2007AGrowth}', '\\cite{Hu2010OnData}', '\\cite{Williams2011TET1Fidelity}'),
    `Expression Data` = c('\\cite{Wang2007AGrowth}', '\\cite{Carroll2006Genome-wideSites}', '\\cite{Williams2011TET1Fidelity}')
) %>%
    xtable(align = 'cllllcc') %>%
    print(floating = FALSE,
          include.rownames = FALSE,
          booktabs = TRUE,
          sanitize.text.function = identity,
          comment = FALSE,
          file = 'article/tables/datasets.tex')

# load required libraries ----
library(tidyverse)

# create directories ----
if(!dir.exists('data/')) dir.create('data/')
if(!dir.exists('data/beta_input/')) dir.create('data/beta_input/')

# download processed data files ----
## reference genomes
c('hg19.refseq',
  'hg19_CTCF_bound.bed',
  'hg18.refseq', 
  'mm9.refseq', 
  'mm9_CTCF_bound.bed',
  'mm10.refseq') %>%
    map(~download.file(
        url = paste0('https://raw.githubusercontent.com/suwangbio/BETA/master/BETA_1.0.7/BETA/references/', .x),
        destfile = paste0('data/beta_input/', .x)
    ))


## bed files
c('3656_peaks.bed',
  '5795_peaks.bed',
  '349_peaks.bed'
) %>%
map(~download.file(
        url = paste0('https://raw.githubusercontent.com/suwangbio/BETA/master/BETA_test_data/', .x),
        destfile = paste0('data/beta_input/', .x)
    ))

## xls files
c('AR_diff_expr.xls',
  'Tet1_diff_expr.xls',
  'ESR1_diff_expr.xls'
) %>%
map(~download.file(
        url = paste0('https://raw.githubusercontent.com/suwangbio/BETA/master/BETA_test_data/', .x),
        destfile = paste0('data/beta_input/', .x)
    ))

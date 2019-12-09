# load required libraries
library(tidyverse)
library(cowplot)

source('code/function-def.R')

# make directories
if(!dir.exists('article/')) dir.create('article/')
if(!dir.exists('article/figures/')) dir.create('article/figures/')

# make figure
plot_grid(plotlist = list(
    plot_target(res_file = 'data/target_output/basic_AR_3656.rds',
                cut_quantiles = c(0, .1, .9, 1)) +
        ggtitle('AR in LNCaP'),
    
    plot_target(res_file = 'data/target_output/basic_ESR1_349.rds',
                cut_quantiles = c(0, .1, .9, 1)) +
        ggtitle('ESR1 in MCF-7'),
    
    plot_target(res_file = 'data/target_output/basic_Tet1_5795.rds',
                cut_quantiles = c(0, .1, .9, 1)) +
        ggtitle('Tet1 in ES cells')
    ),
    nrow = 1,
    scale = .9,
    labels = 'AUTO',
    label_fontface = 'plain',
    label_size = 10) %>%
    ggsave(filename = 'article/figures/ecdf.png',
           width = 24, height = 8, units = 'cm')

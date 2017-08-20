# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# this script contains ggplot theme and color specifications 
# used in plotting figures, as well as some useful functions for plotting


# ----------------------------------------------------------------------
# Loading packages
# ----------------------------------------------------------------------

# packages I usually use for creating figures 
packages <- c('ggplot2', 'dplyr', 'reshape2', 'ggrepel',
              'doParallel', 'foreach', 'doRNG')
lapply(packages, library, character.only = TRUE)


# ----------------------------------------------------------------------
# Fonts and ggplot themes
# ----------------------------------------------------------------------

# font setup
# loadfonts()
fontSetup <- "Helvetica"
fontSize <- 3.15 
pointSize <- 1.5
themeFontSize <- 12

# theme with font sizes adjusted for plotting figures with tikz device
pdftheme <- 
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(lineend = 4, linetype = 1),
        axis.ticks.y = element_line(lineend = 4, linetype = 1),
        axis.ticks = element_line (colour = "black", size = 0.3), 
        axis.text = element_text(size = themeFontSize, colour = "black"),
        axis.text.x = element_text(vjust = 0.5),
        axis.title = element_text(size = themeFontSize + 1),
        axis.title.y = element_text(vjust = 1.8),
        axis.title.x = element_text(vjust = -.8),
        legend.title = element_blank(),
        legend.justification = c(1,0),
        legend.position = c(1,0),
        legend.text = element_text(size = themeFontSize),
        legend.key = element_rect(fill = "#FFFFFF"),
        legend.key.height = unit(0.8,"line"),
        strip.text = element_text(size = themeFontSize + 1),
        strip.background = element_rect(fill = "#FFFFFF"),
        text = element_text(family = fontSetup),
        validate = TRUE
    )

# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# this script contains ggplot theme and color specifications 
# used in plotting figures, as well as some useful functions for plotting


# ----------------------------------------------------------------------
# Loading packages
# ----------------------------------------------------------------------

# packages I usually use for creating figures 
packages <- c('ggplot2', 'dplyr', 'reshape2', 
              'doParallel', 'foreach', 'doRNG')
lapply(packages, library, character.only = TRUE)


# ----------------------------------------------------------------------
# Ggplot themes
# ----------------------------------------------------------------------

# theme with font sizes adjusted for plotting figures with tikz device
pdftheme <- 
    theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.margin = unit(1.3, "lines"),
        axis.ticks.x = element_line(linetype = 0),
        axis.ticks.y = element_line(lineend=4, linetype = 1, colour = "black"),
        axis.ticks = element_line (colour = "black", size = 0.3), 
        axis.text = element_text(size=12, colour = "black"),
        axis.text.x = element_text(vjust = 0.5, size=12),
        axis.title=element_text(size=13),
        axis.title.y = element_text(vjust = 1.8),
        axis.title.x = element_text(vjust = -.8),
        legend.title = element_blank(),
        legend.justification = c(1,0.5),
        legend.position = c(1,0.5),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="#FFFFFF"),
        strip.text = element_text(size=13),
        strip.background = element_rect(fill="#FFFFFF")
    )

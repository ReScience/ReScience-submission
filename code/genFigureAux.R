# ----------------------------------------------------------------------
# Description
# ----------------------------------------------------------------------

# Script for Hinton, Nowlan 1987 article
# It produces additional figure with path for undecided allele for each 
# simulation run, and histogram of end points, with original results
# from the reference paper overlaid 


# ----------------------------------------------------------------------
# Loading the data
# ----------------------------------------------------------------------

# you will need to set the directory if using the script interactively, 
# below you see example of my path to the folder
# setwd("/home/hstojic/Research/replications/las_HintonNowlan1987_repl/code")

# house keeping
rm(list = ls())

# load libraries and auxiliary functions
source("utils.R")

# load the dataset
evoRes <- read.csv(file = "../data/simulation.csv", stringsAsFactors = FALSE)


# ----------------------------------------------------------------------
# Reshaping the data
# ----------------------------------------------------------------------

# reshaping data for plotting
colnames(evoRes) <- c("Simulation", "Correct", "Incorrect", "Undecided")

# number of epochs we plot
noEpochs <- 1000

# data for evolution paths of all simulations for undecided alleles
plotFrame <- evoRes %>% 
    group_by(Simulation) %>% 
    mutate(Epoch = 1:n()) %>%
    mutate(lineLabel = ifelse(Epoch == 600 & Simulation == 2, "Single simulation", NA) )

plotFrame_mean <- plotFrame %>% 
    group_by(Epoch) %>% 
    summarize(Undecided = mean(Undecided)) %>%
    mutate(lineLabel = ifelse(Epoch == 400, "Mean", NA) )

# data for the histogram of the end points
plotHist <- filter(plotFrame, Epoch == max(Epoch))

# the endpoint of the original Hinton & Nowlan allele frequencies 
# aaproximate, from figure 2
undecidedRefArticle <- 0.46
refArticle <- 
    data.frame(Epoch = 1:noEpochs, Undecided = undecidedRefArticle) %>%
    mutate(lineLabel = ifelse(Epoch == 800, "Hinton & nowlan (1987)", NA) )


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

# path figure
relFrequencies_allSims <- 
    ggplot(
        plotFrame, 
        aes(x = Epoch, y = Undecided)) +
    geom_line(aes(group = Simulation), alpha = 0.3) + 
    geom_label_repel( 
        aes(x = Epoch, y = Undecided, label = lineLabel),
        colour = "black", segment.alpha = 1,
        segment.colour = "black", segment.size = 0.5, force = 1,
        min.segment.length = unit(0.1, "lines"),
        size = fontSize,
        inherit.aes = FALSE) + 
    geom_line(
        data = plotFrame_mean, 
        aes(x = Epoch, y = Undecided), color = "red") + 
    geom_label_repel( 
        data = plotFrame_mean, 
        aes(x = Epoch, y = Undecided, label = lineLabel),
        colour = "black", segment.alpha = 1,
        segment.colour = "black", segment.size = 0.5, force = 1,
        min.segment.length = unit(0.1, "lines"),
        size = fontSize,
        inherit.aes = FALSE) + 
    geom_line(data = refArticle, 
        aes(x = Epoch, y = Undecided), color = "blue") +
    geom_label_repel( 
        data = refArticle, 
        aes(x = Epoch, y = Undecided, label = lineLabel),
        colour = "black", segment.alpha = 1,
        segment.colour = "black", segment.size = 0.5, force = 1,
        min.segment.length = unit(0.1, "lines"),
        size = fontSize,
        inherit.aes = FALSE) + 
    scale_x_continuous("\nGenerations",
                     limits = c(1, noEpochs),
                     breaks = c(1, seq(0, noEpochs, 100)[2:length(seq(0, noEpochs, 100))]),
                     labels = c(1, seq(0, noEpochs, 100)[2:length(seq(0, noEpochs, 100))])) +
    scale_y_continuous("Relative frequency of allele\n",
                       limits = c(0, 1),
                       breaks = seq(0, 1, 0.1)) +
    pdftheme + 
    theme(legend.position = "none")


# histogram
endPointHistogram <- 
    ggplot(plotHist, aes(x = Undecided)) +
    geom_histogram(binwidth = 0.01, fill = "white", color = "black") +
    # geom_bar(fill = "white", color = "black") +
    scale_x_continuous("\nFinal relative frequency of Undecided alleles") +
    scale_y_continuous("Count\n",
                       limits = c(0, 30),
                       breaks = seq(0, 30, 5)) +
    pdftheme 


# ----------------------------------------------------------------------
# Saving figures
# ----------------------------------------------------------------------

### PDFs

cairo_pdf("../article/Figure_relFrequencies_allSims.pdf", height = 4, width = 7, onefile = TRUE)
print(relFrequencies_allSims)
dev.off()

cairo_pdf("../article/Figure_endPointHistogram.pdf", height = 4, width = 7, onefile = TRUE)
print(endPointHistogram)
dev.off()


### SVGs

svg("../article/Figure_relFrequencies_allSims.svg", height = 4, width = 7, onefile = TRUE)
print(relFrequencies_allSims)
dev.off()

svg("../article/Figure_endPointHistogram.svg", height = 4, width = 7, onefile = TRUE)
print(endPointHistogram)
dev.off()

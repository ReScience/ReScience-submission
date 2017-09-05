# ----------------------------------------------------------------------
# Description
# ----------------------------------------------------------------------

# Script for Hinton, Nowlan 1987 article
# It produces the figure 2 from the article (note that the caption of
# figure 2 is switched with caption of figure 1)


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

# adding nice names for the data
colnames(evoRes) <- c("Simulation", "Correct", "Incorrect", "Undecided")

# detect the number of simulations in the data
noSim <- max(evoRes$Simulation)

# computing means
plotFrame_means <- evoRes %>% 
    group_by(Simulation) %>% 
    mutate(Epoch = 1:n()) %>%
    group_by(Epoch) %>% 
    summarize(
        Correct = mean(Correct),
        Incorrect = mean(Incorrect),
        Undecided = mean(Undecided))  %>%
    melt( 
        id.vars = "Epoch",
        measure.vars = c("Correct", "Incorrect", "Undecided"),
        value.name = "proportion_mean",
        variable.name = "alleleType")

# computing standard errors
plotFrame_ses <- evoRes %>% 
    group_by(Simulation) %>% 
    mutate(Epoch = 1:n()) %>%
    group_by(Epoch) %>% 
    summarize(
        Correct = sd(Correct)/sqrt(noSim),
        Incorrect = sd(Incorrect)/sqrt(noSim),
        Undecided = sd(Undecided)/sqrt(noSim)) %>%
    melt( 
        id.vars = "Epoch",
        measure.vars = c("Correct", "Incorrect", "Undecided"),
        value.name = "proportion_se",
        variable.name = "alleleType") 

plotFrame <- left_join(plotFrame_means, plotFrame_ses, 
    by = c("Epoch", "alleleType")) 


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

# plot of all 1000 epochs
noEpochs <- 1000
relFrequencies1000 <- 
    ggplot(
        mutate(plotFrame, 
            lineLabel = ifelse(Epoch == 0.8*noEpochs, 
                as.character(alleleType), NA)), 
        aes(x = Epoch, y = proportion_mean)) +
    geom_line(aes(group = alleleType, linetype = alleleType)) + 
    geom_label_repel( 
        aes(x = Epoch, y = proportion_mean, label = lineLabel),
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
    scale_linetype_manual("",
                          values = c(3,2,1),
                          labels = c("Correct alleles", 
                                     "Incorrect alleles",
                                     "Undecided alleles")) +
    pdftheme + 
    theme(legend.position = "none")

# if there is more than 1 simulation we plot confidence interval as well
if (noSim > 1) {
    relFrequencies1000 <- relFrequencies1000 +  
        geom_ribbon(
            aes(ymin = ifelse((proportion_mean - proportion_se) < 0, 0,
                              proportion_mean - proportion_se), 
                ymax = ifelse((proportion_mean + proportion_se) > 1, 1, 
                              proportion_mean + proportion_se), 
            group = alleleType), alpha = 0.2) 
}


# plot of first 50 epochs, same as in the article
noEpochs <- 50
relFrequencies50 <- 
    ggplot(
        filter(plotFrame, Epoch <= noEpochs) %>%
        mutate(lineLabel = ifelse(Epoch == 0.8*noEpochs, 
            as.character(alleleType), NA)), 
        aes(x = Epoch, y = proportion_mean)) +
    geom_line(aes(group = alleleType, linetype = alleleType)) + 
    geom_label_repel( 
        aes(x = Epoch, y = proportion_mean, label = lineLabel),
        colour = "black", segment.alpha = 1,
        segment.colour = "black", segment.size = 0.5, force = 1,
        min.segment.length = unit(0.1, "lines"),
        size = fontSize,
        inherit.aes = FALSE) + 
    scale_x_continuous("\nGenerations",
                     limits = c(1, noEpochs),
                     breaks = c(1, seq(0, noEpochs, 10)[2:length(seq(0, noEpochs, 10))]),
                     labels = c(1, seq(0, noEpochs, 10)[2:length(seq(0, noEpochs, 10))])) +
    scale_y_continuous("Relative frequency of allele\n",
                       limits = c(0, 1),
                       breaks = seq(0, 1, 0.1)) +
    scale_linetype_manual("",
                          values = c(3,2,1),
                          labels = c("Correct alleles", 
                                     "Incorrect alleles",
                                     "Undecided alleles")) +
    pdftheme + 
    theme(legend.position = "none")

# if there is more than 1 simulation we plot confidence interval as well
if (noSim > 1) {
    relFrequencies50 <- relFrequencies50 +  
        geom_ribbon(
            aes(ymin = ifelse((proportion_mean - proportion_se) < 0, 0,
                              proportion_mean - proportion_se), 
                ymax = ifelse((proportion_mean + proportion_se) > 1, 1, 
                              proportion_mean + proportion_se), 
            group = alleleType), alpha = 0.2) 
}


# ----------------------------------------------------------------------
# Saving figures
# ----------------------------------------------------------------------

### PDFs

cairo_pdf("../article/Figure2.pdf", height = 4, width = 7, onefile = TRUE)
print(relFrequencies50)
dev.off()

cairo_pdf("../article/Figure2_1000.pdf", height = 4, width = 7, onefile = TRUE)
print(relFrequencies1000)
dev.off()


### SVGs

svg("../article/Figure2.svg", height = 4, width = 7, onefile = TRUE)
print(relFrequencies50)
dev.off()

svg("../article/Figure2_1000.svg", height = 4, width = 7, onefile = TRUE)
print(relFrequencies1000)
dev.off()
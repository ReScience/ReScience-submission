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

# load libraries and auxiliary functions
source("utils.R")

# load the dataset
evoRes <- read.csv(file = "../data/simulation.csv", stringsAsFactors = FALSE)


# ----------------------------------------------------------------------
# Reshaping the data
# ----------------------------------------------------------------------

# reshaping data for plotting
colnames(evoRes) <- c("Simulation", "Correct", "Incorrect", "Undecided")

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

plotFrame_ses <- evoRes %>% 
    group_by(Simulation) %>% 
    mutate(Epoch = 1:n()) %>%
    group_by(Epoch) %>% 
    summarize(
        Correct = sd(Correct)/100,
        Incorrect = sd(Incorrect)/100,
        Undecided = sd(Undecided)/100) %>%
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
            lineLabel = ifelse(Epoch == 700, as.character(alleleType), NA)), 
        aes(x = Epoch, y = proportion_mean)) +
    geom_line(aes(group = alleleType, linetype = alleleType)) + 
    geom_ribbon(
        aes(ymin = ifelse((proportion_mean - proportion_se) < 0, 0,
                          proportion_mean - proportion_se), 
            ymax = ifelse((proportion_mean + proportion_se) > 1, 1, 
                          proportion_mean + proportion_se), 
        group = alleleType)) + 
    geom_label_repel( 
        aes(x = Epoch, y = proportion_mean, label = lineLabel),
        colour = "black", 
        segment.color = "black", segment.size = 0.2, force = 10,
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


# plot of first 50 epochs, same as in the article
noEpochs <- 50
relFrequencies50 <- 
    ggplot(
        filter(plotFrame, Epoch <= noEpochs) %>%
        mutate(lineLabel = ifelse(Epoch == 45, as.character(alleleType), NA)), 
        aes(x = Epoch, y = proportion_mean)) +
    geom_line(aes(group = alleleType, linetype = alleleType)) + 
    geom_ribbon(
        aes(ymin = ifelse((proportion_mean - proportion_se) < 0, 0,
                          proportion_mean - proportion_se), 
            ymax = ifelse((proportion_mean + proportion_se) > 1, 1, 
                          proportion_mean + proportion_se), 
        group = alleleType)) + 
    geom_label_repel( 
        aes(x = Epoch, y = proportion_mean, label = lineLabel),
        colour = "black", 
        segment.color = "black", segment.size = 0.2, force = 10,
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


# ----------------------------------------------------------------------
# Saving figures
# ----------------------------------------------------------------------

cairo_pdf("../article/Figure2.pdf", height = 4, width = 7, onefile = TRUE)
print(relFrequencies50)
dev.off()

cairo_pdf("../article/Figure2_1000.pdf", height = 4, width = 7, onefile = TRUE)
print(relFrequencies1000)
dev.off()

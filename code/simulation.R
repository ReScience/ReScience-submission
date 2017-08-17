# ----------------------------------------------------------------------
# Description
# ----------------------------------------------------------------------

# Script for running the simulation from Hinton, Nowlan 1987 article


# ----------------------------------------------------------------------
# Setting the parameters
# ----------------------------------------------------------------------

# you will need to set the directory if using the script interactively, 
# below you see example of my path to the folder
# setwd("/home/hstojic/Research/replications/las_HintonNowlan1987_repl/code")

# load libraries and auxiliary functions
source("utils.R")

# we set base probabilities  
# for elements in the genome, according to the article
baseProbs <- c(0.25, 0.25, 0.5)  # 0, 1, NA (? from the article)
noGenerations <- 1000
noAgents <- 1000
lifetime <- 1000
noSim <- 100
seed <- 1234

# to speed up simulations, we run it in parallel 
noCores <-  ifelse(detectCores() == 1, 1, detectCores() - 1)
registerDoParallel(noCores)
registerDoRNG(seed)  


# ----------------------------------------------------------------------
# Defining functions
# ----------------------------------------------------------------------

# function for an agent going through all the rounds in its lifetime
# trying to learn the target
agent <- function(genotype, target, lifetime = 1000) {

    # identify how much evolution left to an agent to learn
    learningSpace <- is.na(genotype)
    learningGoal <- target[learningSpace]

    # if evolution did not leave anything for an agent to learn
    # we return back the original genotype and lifetime 
    if (!any(learningSpace)) {
        return(list(trialSolved = lifetime, 
                    genome = genotype))
    
    # if evolutionary learning is incorrect, no amount of individual 
    # learning will yield any benefit, we end the lifetime
    } else if (!all(genotype[!learningSpace] == target[!learningSpace])) {
        learnSolution <- learn(learningGoal,1)$solution
        genotype[learningSpace] <- learnSolution
        return(list(trialSolved = lifetime, 
                    genome = genotype))
    
    # else agent does individual learning
    } else {
        learningOutcome <- learn(learningGoal, lifetime)
        genotype[learningSpace] <- learningOutcome$solution
        return(list(trialSolved = learningOutcome$trialSolved, 
                    genome = genotype))
    }
}

# auxiliary function for performing the learning
learn <- function(learningGoal, lifetime) {
    
    # generate candidate solution at random
    learningAmount <- length(learningGoal)
    candidates <- matrix(rbinom(learningAmount*lifetime, 1, 0.5), 
                         nr = learningAmount)
    candidatesEvaluated <- candidates == learningGoal

    # detect trial in which solution appeared, if at all
    solutions <- apply(candidatesEvaluated, 2, all)
    if(any(solutions)) {
        trialSolved <- min(which(solutions)) 
    } else {
        trialSolved <- lifetime
    }
    solution <- candidates[,trialSolved]
    return(list(trialSolved = trialSolved, solution = solution))
}

# function for producing a child genome from two parent genomes
genChild <- function(genome1, genome2) {
    N <- length(genome1)
    crossover <-sample((0:N) + 0.5, 1)
    genomeChild <- c(genome1[1:N < crossover], genome2[1:N > crossover])
    return(genomeChild)
}

# function for evaluating the fitness of a genome
# here depends on how fast it learns the solution in its lifetime
evaluateFitness <- function(lifetime, trialSolved) {
    1 + (19 * (lifetime - trialSolved))/1000
}

# main function that performs evolution of the population of agents
evolveAgents <- function(target, baseProbs, 
                         noAgents = 1000, lifetime = 1000, 
                         noGenerations = 1000, 
                         diagnostics = TRUE, seed) {
    
    # set the seed for reproducability
    if (missing(seed)) {
        seed <- round(runif(1)*10000)
    }
    set.seed(seed)

    # initialize  agents' genomes and empty frame for results
    genomeLength <- length(target)
    genomes <- matrix(sample(c(0,1,NA), genomeLength * noAgents, 
                             prob = baseProbs, replace = TRUE),
                      nc = noAgents)
    alleleFreqMean <- matrix(NA, nr = noGenerations, nc = 3)

    for (gen in 1:noGenerations) {
        if (diagnostics) print(gen)
        
        genomeFitness <- rep(NA, noAgents)
        alleleFreq <- matrix(NA, nr = noAgents, nc = 3)


        ### learning
        
        for (ag in 1:noAgents) {
            genotype <- genomes[,ag]
            learningOutcome <- agent(genotype, target, lifetime)
            
            # compute mating score
            genomeFitness[ag] <- 
                evaluateFitness(lifetime, learningOutcome$trialSolved)

            # frequency of correct, incorrect and undecided alleles
            correct <-  genotype == target
            alleleFreq[ag,] <- c(
                sum(correct, na.rm = TRUE), 
                sum(!correct, na.rm = TRUE), 
                sum(is.na(correct))
                ) / genomeLength
        }

        # final result of the whole generation, mean allele frequencies across
        # all genomes
        alleleFreqMean[gen,] <- colMeans(alleleFreq)


        ### mating
        
        # first we compute the probabilities from genomeFitnesss
        # then we do as many matings as there are agents
        # for each mating parent is unique
        # with genChild function we produce a child genome
        reproductionProb <- genomeFitness / sum(genomeFitness)
        newGenomes <- matrix(NA, nr = genomeLength, nc = noAgents)
        for (mating in 1:noAgents) {
            pair <- sample(1:noAgents, 2, prob = reproductionProb, 
                replace = FALSE) 
            newGenomes[,mating] <- genChild(genomes[,pair[1]], 
                                            genomes[,pair[2]])
        } 

        # finally, update the genomes for the next generation
        genomes <- newGenomes
    }

    return(alleleFreqMean)
}



# ----------------------------------------------------------------------
# Simulation
# ----------------------------------------------------------------------

# evolving population of agents, repeating the simulation for noSim times
results <- foreach (sim = 1:noSim, .combine = rbind) %dorng% {
    print(sim)

    # we pick a targeted pattern of alleles - pattern which agents need
    # to learn, any will do 
    target <- rbinom(20,1,runif(1))
    print(target)

    # single simulation of an evolution
    evoRes <- evolveAgents(target, baseProbs, 
                           noAgents = noAgents, lifetime = lifetime, 
                           noGenerations = noGenerations, 
                           diagnostics = FALSE)
    return(data.frame(sim, evoRes))
}

# saving the data
write.csv(results, 
          file = "../data/simulation.csv",
          row.names = FALSE)


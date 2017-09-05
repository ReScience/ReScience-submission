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

# house keeping
rm(list = ls())

# load libraries and auxiliary functions
source("utils.R")

# we set base probabilities  
# for elements in the genome, according to the article
noAlleles <- 20
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
# trying to learn the targeted genome, gstar
# it is used in individual learning part of the loop from Algorithm 1
agentLearning <- function(g, gstar, lifetime = 1000) {

    # identify part of the genome, g, that evolution left to an agent to learn
    glearnID <- is.na(g)
    gstarL <- gstar[glearnID]

    # if evolution did not leave anything for an agent to learn
    # we return back the original genotype and lifetime 
    if (!any(glearnID)) {
        return(list(stepSolved = lifetime, genome = g))
    
    # if evolutionary learning is incorrect, no amount of individual 
    # learning will yield any benefit, we end the lifetime
    } else if (!all(g[!glearnID] == gstar[!glearnID])) {
        learnSolution <- learning(gstarL,1)$solution
        g[glearnID] <- learnSolution
        return(list(stepSolved = lifetime, genome = g))
    
    # else agent does individual learning
    } else {
        learningOutcome <- learning(gstarL, lifetime)
        g[glearnID] <- learningOutcome$solution
        return(list(stepSolved = learningOutcome$stepSolved, 
                    genome = g))
    }
}

# auxiliary function for performing the individual learning 
# for a single agent
# it is used in individual learning part of the loop from Algorithm 1
# (step = 1 : lifetime)
learning <- function(gstarL, lifetime) {
    
    # generate candidate solution randomly
    candidates <- matrix(
        rbinom(length(gstarL)*lifetime, 1, 0.5), 
        nr = length(gstarL)
    )
    candidatesEvaluated <- candidates == gstarL

    # detect step in which solution appeared, if at all
    solutions <- apply(candidatesEvaluated, 2, all)
    if(any(solutions)) {
        stepSolved <- min(which(solutions)) 
    } else {
        stepSolved <- lifetime
    }
    solution <- candidates[,stepSolved]
    return(list(stepSolved = stepSolved, solution = solution))
}

# function for producing a child genome from two parent genomes
# used in the loop from Algorithm 1 that generates children genomes 
genChild <- function(genome1, genome2) {
    N <- length(genome1)
    crossover <-sample((0:N) + 0.5, 1)
    genomeChild <- c(genome1[1:N < crossover], genome2[1:N > crossover])
    return(genomeChild)
}

# function for evaluating the fitness of a genome
# fitness depends on how fast it learns the solution in its lifetime
# used in the loop from Algorithm 1 that does individual learning 
evaluateFitness <- function(lifetime, stepSolved, noAlleles) {
    1 + ((noAlleles - 1) * (lifetime - stepSolved))/lifetime
}

# main function that performs evolution of the population of agents
# performs the second loop from Algorithm 1 (gen = 1 : noGenerations)
evolveAgents <- function(noAlleles, gstar, baseProbs, 
                         noAgents = 1000, lifetime = 1000, 
                         noGenerations = 1000, 
                         diagnostics = TRUE, seed) {
    
    # set the seed for reproducability
    if (missing(seed)) {
        seed <- round(runif(1)*10000)
    }
    set.seed(seed)

    # initialize  agents' genomes and empty frame for results
    gLength <- length(gstar)
    G <- matrix(
        sample(c(0,1,NA), gLength*noAgents, prob = baseProbs, replace = TRUE),
        nc = noAgents
    )
    alleleFreqMean <- matrix(NA, nr = noGenerations, nc = 3)

    for (gen in 1:noGenerations) {
        if (diagnostics) print(gen)
        
        f <- rep(NA, noAgents)
        alleleFreq <- matrix(NA, nr = noAgents, nc = 3)


        ### individual learning
        
        # this is the loop from Algorithm 1 that does individual learning
        for (ag in 1:noAgents) {
            g <- G[,ag]
            learningOut <- agentLearning(g, gstar, lifetime)
            
            # compute fitness scores
            f[ag] <- 
                evaluateFitness(lifetime, learningOut$stepSolved, noAlleles)

            # frequency of correct, incorrect and undecided alleles
            correct <-  g == gstar
            alleleFreq[ag,] <- c(
                sum(correct, na.rm = TRUE), 
                sum(!correct, na.rm = TRUE), 
                sum(is.na(correct))
                ) / gLength
        }

        # final result of the whole generation, mean allele frequencies across
        # all genomes
        alleleFreqMean[gen,] <- colMeans(alleleFreq)


        ### evolutionary learning
        
        # first we compute the probabilities from f
        # then we do as many matings as there are agents
        # for each mating parent is unique
        # with genChild function we produce a child genome
        p <- f/sum(f)
        Gnew <- matrix(NA, nr = gLength, nc = noAgents)

        # this is the loop from Algorithm 1 that generates children genomes 
        for (i in 1:noAgents) {
            pair <- sample(1:noAgents, 2, prob = p, replace = FALSE) 
            Gnew[,i] <- genChild(G[,pair[1]], G[,pair[2]])
        } 

        # finally, update the genomes for the next generation
        G <- Gnew
    }

    return(alleleFreqMean)
}



# ----------------------------------------------------------------------
# Simulation
# ----------------------------------------------------------------------

# evolving population of agents, repeating the simulation for noSim times
# this is the first loop from Algorithm 1 
results <- foreach (sim = 1:noSim, .combine = rbind) %dorng% {

    # we pick a targeted pattern of alleles - pattern which agents need
    # to learn, any will do 
    gstar <- rbinom(noAlleles, 1, runif(1))

    # single simulation of an evolution
    startTime <-  Sys.time()
    evoRes <- evolveAgents(noAlleles, gstar, baseProbs, 
                           noAgents = noAgents, lifetime = lifetime, 
                           noGenerations = noGenerations, 
                           diagnostics = FALSE)
    endTime <-  Sys.time()
    runTime <- endTime - startTime
    cat("Simulation: ", sim, " | Execution time: ", runTime, "min\n")

    # provide some estimate of remaining runs duration
    cat("Estimated remaining time: ", 
        as.numeric(runTime)*(noSim-sim)/noCores, "min\n")

    return(data.frame(sim, evoRes))
}

# saving the data
write.csv(results, 
          file = "../data/simulation.csv",
          row.names = FALSE)


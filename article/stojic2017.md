---
Title: "How learning can guide evolution"
Author:
  - name: Hrvoje Stojić
    affiliation: 1
Address:
  - code:    1
    address: Department of Economics and Business, Universitat Pompeu Fabra, Barcelona, Spain
Contact:
  - hrvoje.stojic@protonmail.com
Editor:
  - Name Surname
Reviewer:
  - Name Surname
  - Name Surname
Publication:
  received:  Aug, 20, 2017
  accepted:  Sep, 1, 2015
  published: Sep, 1, 2015
  volume:    "**1**"
  issue:     "**1**"
  date:      Sep 2015
Repository:
  article:   "https://github.com/ReScience-Archives/Stojic-2017/tree/master/article"
  code:      "https://github.com/ReScience-Archives/Stojic-2017/tree/master/code"
  data:      "https://github.com/ReScience-Archives/Stojic-2017/tree/master/data"
  notebook:  
Reproduction:
  - "*How learning can guide evolution*, G. E. Hinton, and S. J. Nowlan, Complex Systems, 1 (3), 495-502, 1987."
Bibliography:
  stojic2017.bib

---

# Introduction

The Lamarckian hypothesis that adaptations accumulated during the lifetime of an individual organism are directly transmitted to the next generation is considered to be incorrect. The reference article [@hinton1987learning] argues that this does not necessarily mean that individual learning cannot exert influence on evolution in an indirect fashion, providing a simple but convincing computational example of such an interaction. The authors use an extreme scenario where the optimization landscape is completely flat except for a single peak. In such an environment, curvature of the optimization surface does not provide any guidance and evolutionary search alone cannot find the maximum. Through a simulation the authors show that such a hostile environment can be tackled by a combination of evolutionary and individual learning. In addition, if organisms can search during their lifetime, those organisms whose genomes are closer to the targeted peak are going to be able to find the peak with individual search. Such organisms will have higher fitness and transmit their genes to the next generation. In effect, capacity for individual learning "creates a hill" to the peak that evolutionary search can climb, as illustrated in Figure 1 in @hinton1987learning^[Note that captions of the figures are switched, caption of Figure 1 is incorrectly placed under Figure 2, and vice versa]. This is one of the clearest illustrations of the benefits of individual learning for evolution and the interaction between the two learning processes.^[An indirect influence of individual learning on evolutionary learning is often called Baldwin effect.] This result has made a substantial impact in the cognitive science literature.

The main aim of this article is to replicate the simulations reported in the reference article [@hinton1987learning]. I have corresponded with one of the authors and learned that the original implementation of the simulations is no longer available. I am also not aware of other implementations elsewhere. The simulations are relatively simple, and I thus propose a replication based on description from the reference article. Replication code relies on R programming language [@R] and several R packages [@ggplot2; @dplyr; @reshape2; @doParallel; @foreach; @doRNG; @ggrepel].


# Methods

A description of the simulations can be found in the caption of Figure 1 in the reference article [@hinton1987learning]. I consider the description to be detailed enough for the implementation. I have followed it closely, with two exceptions. First, I have repeated the evolutionary search many times to average out potential noise. In the reference article, such iterations are not explicitly mentioned, and moreover, when interpreting their results from Figure 2, the authors state that Figure 2 presents results of "... a typical evolutionary search ...", suggesting it is a single run of the simulation. Second, I let the population of agents evolve for 1000 generations instead of 50 generations reported in @hinton1987learning. Reason for this will become clear in the Results section.

I describe my implementation in detail in Algorithm 1. Parameters for the simulation are summarized in Table @tbl:parameters.

\begin{algorithm}
\caption{Simulation description}
    \begin{algorithmic}[1]
    \For{$sim=1:noSim$}
        \State{Generate targeted genome, $g^*=\{a_i\}_{i=1}^{noAlleles}$: $a_i \sim Bernoulli(p)$, }
        \State {$p \sim Uniform(0,1)$ }
        \State{Generate initial genomes of the agents, $G=\{g_i\}_{i=1}^{noAgents}$: $\forall g_i$ generate} 
        \State {$\{a_i\}_{i=1}^{noAlleles}$ by sampling with replacement from $\{0,1,?\}$ according to $p_0,p_1,p_?$.}
        \For{$gen=1:noGenerations$}
            \For{$ag=1:noAgents$} \Comment {\textbf{Individual learning}}
                \For{$step=1:lifetime$}
                    \State {For all ? alleles in agent's genome: $a_i \sim Bernoulli(0.5)$}
                    \If{$g_i == g^*$ after individual learning} \State {\textbf{break}} \EndIf
                \EndFor
                \State {Evaluate fitness: $f_{ag} = 1 - 19(lifetime-step)/1000$}
                \State {Record frequency of Correct, Incorrect and Undecided alleles.}
            \EndFor
            \State Compute parenting probabilities: 
            \State $p_{ag}=f_{ag}/\sum_{i=1}^{noAgents}f_{i}$, $\forall ag$
            \For{$i=1:|G|$} \Comment {\textbf{Generate children genomes}}
                \State {Choose two parents by sampling with replacement according to $p_{ag}$}
                \State {Generate new genome, $g_i^{child}$: randomly choose a cross-over point, 
                \State copy all alleles from the first parent up to the cross-over point, 
                \State and from the second parent beyond the cross-over point.}
            \EndFor
            \State {Update genomes: $g_i \gets g_i^{child}$, $\forall i \in G$}
        \EndFor
    \EndFor
    \end{algorithmic}
\end{algorithm}



Table                Value      Description
-------------------- ---------- -------------------------------------------
noSim                100        Number of times evolutionary search is repeated  
noGenerations        1000       Number of generations for evolutionary  algorithm  
noAgents             1000       Number of agents in a population  
lifetime             1000       Number of cycles available for individual learning  
noAlleles            20         Number of alleles in agent's genome  
$p_0$                0.25       Probability that allele is of type zero  
$p_1$                0.25       Probability that allele is of type one  
$p_?$                0.50       Probability that allele is of type ?
-------------------- ---------- -------------------------------------------
Table: Simulation parameters. {#tbl:parameters}


# Results

The simulation results are summarized in Figure @fig:relFrequencies50, which corresponds to Figure 2 from the reference article. Qualitatively, the results are very similar - there is an obvious increase in the proportion of alleles correctly set by the evolutionary learning (close to 75% by generation 50) and proportion of incorrect alleles decreases to zero by generation 50. This is clear evidence of individual learning guiding evolutionary search, even though there is no direct transmission of knowledge acquired during individual learning from one generation to the next.


![The evolution of the relative frequencies of the three possible
types of allele. Correct alleles are those that are set by evolutionary learning and correspond to the targeted pattern. Incorrect alleles represent the proportion of incorrectly set alleles, while Undecided represent the proportion of alleles left for individual learning. Proportions are means of 100 simulation runs, and barely visible grey ribbons are standard errors. This is the reproduction of Figure 2 in the reference article [@hinton1987learning].](Figure2.pdf){#fig:relFrequencies50}


There is an important difference with respect to the reference article, however. In the original simulations, the proportion of Undecided alleles that are left to individual learning is relatively constant, close to initial 50%. The authors point this out as an interesting result, suggesting that there is little selective pressure to specify last few connections by evolutionary learning because those can be quickly set through individual learning. In contrast, the results in Figure @fig:relFrequencies50 show that the proportion of Undecided alleles is continuously decreasing. In fact, whereas the proportion of Correct alleles in the reference article flattens out, in my case it continues increasing, with the increase coming from the proportion of Undecided alleles decreasing further. To make sure that this is really a long run trend I have let populations evolve for a larger number of generations than in the reference article. Results of all 1000 generations are presented in Figure @fig:relFrequencies1000. By generation 1000 the proportion of Undecided alleles drops to 9%. Even though the decrease is smaller with each generation, it does not seem to stop and flatten out by the end of the simulation.

What is the source of this particular difference in results? One explanation is that the results reported in Figure 2 in the reference article stem from a non-typical draw. This is unlikely, as upon examining the evolutionary paths of Undecided alleles in all simulation runs, they all show a strong downward trend already in the first 50 generations, differing only in the exact generation when this decrease onsets. Another explanation is that the original implementation is flawed as there actually is a selection pressure for eliminating individual learning after some time. The fitness function described in the reference article is such that the smaller the number of steps an agent takes in an individual lifetime to reach the peak, the higher the probability of mating and transmitting genes to the next generation. The number of steps is a direct function of the number of undecided alleles in the genome. Hence, there is a selection pressure for decreasing the number of Undecided alleles. In personal correspondence, one of the authors agreed with this interpretation, stating that indeed there should be a selective pressure for a decrease in number of Undecided alleles. The pressure should decrease as a power function of the mean number of steps an agent needs to reach the peak, which explains rapidly decreasing improvements with each generation.


![The evolution of the relative frequencies of the three possible
types of allele, with a larger number of generations. Proportions are means of 100 simulation runs, and barely visible grey ribbons are standard errors.](Figure2_1000.pdf){#fig:relFrequencies1000}



# Conclusion

Qualitatively, the result is comparable to that reported in the reference article, but there are relatively important differences. The proportion of Undecided alleles in the implementation proposed here diminish over time, whereas in the original implementation they stay close to the initial values. I attribute the difference to a flaw in the original implementation, providing arguments based on a fitness function proposed in the reference article. Although the main result of the reference article - that individual learning can indirectly guide evolutionary search - has been reproduced, the difference in simulation results suggests that individual learning might be eradicated over time. This decrease in usefulness of individual learning occurs in static environments where targeted genome stays fixed, as in these simulations. In more realistic scenarios, where environments are non-stationary and the targeted genome would be changing over time, individual learning would have a more stable role.


# References

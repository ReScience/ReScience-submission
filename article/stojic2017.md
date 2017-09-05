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

The Lamarckian hypothesis that adaptations accumulated during the lifetime of an individual organism are directly transmitted to the next generation is considered to be incorrect. The reference article [@hinton1987learning] argues that this does not necessarily mean that individual learning cannot exert influence on evolution in an indirect fashion, providing a simple but convincing computational example of such an interaction. The authors use an extreme scenario where the optimization landscape is completely flat except for a single peak. In such an environment, curvature of the optimization surface does not provide any guidance and evolutionary search alone cannot find the maximum. Through a simulation the authors show that such a hostile environment can be tackled by a combination of evolutionary and individual learning. For organisms that can search during their lifetime, those organisms whose genomes are closer to the targeted peak are going to be able to find the peak with individual search. Such organisms will have higher fitness and transmit their genes to the next generation. In effect, capacity for individual learning "creates a hill" to the peak that evolutionary search can climb, as illustrated in Figure 1 in @hinton1987learning^[Note that captions of the figures are switched in the reference article, caption of Figure 1 is incorrectly placed under Figure 2, and vice versa.]. In contrast to the Lamarckian hypothesis, this is an indirect transmission of individual learning on evolutionary learning, what is often called the Baldwin effect [@Baldwin1896]. This is one of the clearest illustrations of the benefits of individual learning for evolution and the interaction between the two learning processes. This article is well regarded and widely cited in cognitive science literature.

I replicated the simulations reported in the reference article [@hinton1987learning] in R [@R], using several packages [@devtools; @ggplot2; @dplyr; @reshape2; @doParallel; @foreach; @doRNG; @ggrepel]. In correspondence with Geoffrey Hinton I learned that the original implementation is not available. 


# Methods

A description of the simulations can be found in the caption of Figure 1 in the reference article [@hinton1987learning]. I consider the description to be detailed enough for the implementation. I have followed it closely, with two exceptions. First, I have repeated the evolutionary search many times to average out potential noise. In the reference article, such iterations are not explicitly mentioned, and moreover, when interpreting their results from Figure 2, the authors state that Figure 2 presents results of "... a typical evolutionary search ...", suggesting it is a single run of the simulation. Second, I let the population of agents evolve for 1000 generations instead of 50 generations reported in @hinton1987learning. Reason for this will become clear in the Results section.

Pseudocode of my implementation is shown in Algorithm 1. Parameters for the simulation are summarized in Table @tbl:parameters.

\begin{algorithm}
\caption{Simulation description}
    \begin{algorithmic}[1]
    \For{$sim=1:noSim$}
        \State{Generate targeted genome, $g^*=\{a_i\}_{i=1}^{noAlleles}$: $a_i \sim Bernoulli(p)$, }
        \State {$p \sim Uniform(0,1)$ }
        \State{Generate initial genomes of the agents, $G=\{g_i\}_{i=1}^{noAgents}$: for all $g_i$ generate} 
        \State {$\{a_i\}_{i=1}^{noAlleles}$ by sampling with replacement from $\{0,1,?\}$ according to $p_0,p_1,p_?$.}
        \For{$gen=1:noGenerations$}
            \For{$ag=1:noAgents$} \Comment {\textbf{Individual learning}}
                \For{$step=1:lifetime$}
                    \State {For all ? alleles in agent's genome: $a_i \sim Bernoulli(0.5)$}
                    \If{$g_i == g^*$ after individual learning} \State {\textbf{break}} \EndIf
                \EndFor
                \State {Evaluate fitness: $f_{ag} = 1 + 19(lifetime-step)/1000$}
                \State {Record frequency of Correct, Incorrect and Undecided alleles.}
            \EndFor
            \State Compute parenting probabilities: 
            \State $p_{ag}=f_{ag}/\sum_{i=1}^{noAgents}f_{i}$, for every agent $ag$
            \For{$i=1:noAgents$} \Comment {\textbf{Generate children genomes}}
                \State {Choose two parents by sampling with replacement according to $p_{ag}$}
                \State {Generate new genome, $g_i^{child}$: randomly choose a cross-over point, 
                \State copy all alleles from the first parent up to the cross-over point, 
                \State and from the second parent beyond the cross-over point.}
            \EndFor
            \State {Update genomes: $g_i \gets g_i^{child}$, for every genome in $G$}
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
$p_0$                0.25       Probability that allele's value is 0  
$p_1$                0.25       Probability that allele's value is 1  
$p_?$                0.50       Probability that allele's value is ? (learnable)  
-------------------- ---------- -------------------------------------------
Table: Simulation parameters. {#tbl:parameters}


# Results

The simulation results are summarized in Figure @fig:relFrequencies50, which corresponds to Figure 2 from the reference article. Qualitatively, the results are very similar - there is an obvious increase in the proportion of alleles correctly set by the evolutionary learning (close to 75% by generation 50) and proportion of incorrect alleles decreases to zero by generation 50. This is clear evidence of individual learning guiding evolutionary search, even though there is no direct transmission of knowledge acquired during individual learning from one generation to the next.


![The evolution of the relative frequencies of the three possible
types of allele. Correct alleles are those that are set by evolutionary learning and correspond to the targeted pattern. Incorrect alleles represent the proportion of incorrectly set alleles, while Undecided represent the proportion of alleles left for individual learning. Proportions are means of 100 simulation runs, and grey shaded area around the mean line represent standard errors of the means. This is the reproduction of Figure 2 in the reference article [@hinton1987learning].](Figure2.pdf){#fig:relFrequencies50}


There is an important difference with respect to the reference article, however. In the original simulations, the proportion of Undecided alleles that are left to individual learning is relatively constant, close to initial 50%. The authors point this out as an interesting result, suggesting that there is little selective pressure to specify last few connections by evolutionary learning because those can be quickly set through individual learning. In contrast, the results in Figure @fig:relFrequencies50 show that the proportion of Undecided alleles is continuously decreasing. In fact, whereas the proportion of Correct alleles in the reference article flattens out, in my case it continues increasing, with the increase coming from the proportion of Undecided alleles decreasing further. To make sure that this is really a long run trend I have let populations evolve for a larger number of generations than in the reference article. Results of all 1000 generations are presented in Figure @fig:relFrequencies1000. By generation 1000 the mean proportion of Undecided alleles drops to 9%, with the decrease being smaller with each generation.


![The evolution of the relative frequencies of the three possible
types of allele, with a larger number of generations. Proportions are means of 100 simulation runs, and grey shaded area around the mean line represent standard errors of the means.](Figure2_1000.pdf){#fig:relFrequencies1000}


What is the source of this particular difference in results? The results reported in Figure 2 in the reference article might stem from a non-typical draw. This is unlikely, as illustrated in Figure @fig:Figure_relFrequencies_allSims. Evolutionary paths of Undecided alleles in all simulation runs show a strong downward trend already in the first 50 generations, differing only in the exact generation when this decrease onsets. As shown in Figure @fig:Figure_endPointHistogram, distribution of proportions of Undecided alleles at the end of the simulation (1000 generations) ranges from 0.05 to 0.20, much lower than 0.46 that was obtained in the reference article.


![The evolution of relative frequencies of Undecided alleles for all 100 simulations, red line marks the mean across all simulations, same as in Figure @fig:relFrequencies1000, while blue line marks approximately the end point of relative frequency from the reference article.](Figure_relFrequencies_allSims.pdf){#fig:Figure_relFrequencies_allSims}


![Histogram of relative frequencies of Undecided alleles at the end of simulations (1000 generations).](Figure_endPointHistogram.pdf){#fig:Figure_endPointHistogram}


The decrease in Undecided alleles is actually expected, once the fitness function described in the reference article is examined more closely. The function is such that the smaller the number of steps an agent takes in an individual lifetime to reach the peak, the higher the probability of mating and transmitting genes to the next generation. The number of steps is a direct function of the number of undecided alleles in the genome. In other words, there exists a selection pressure for decreasing the number of Undecided alleles. In personal correspondence, Geoffrey Hinton agreed with this argument, stating that the selective pressure should decrease as a power function of the mean number of steps an agent needs to reach the peak. Computational power at the time was orders of magnitude smaller than today, preventing them to explore the space more thoroughly.


# Conclusion

Qualitatively, the result is comparable to that reported in the reference article, but there are relatively important differences. The proportion of Undecided alleles in the implementation proposed here diminish over time, whereas in the original implementation they stay close to the initial values. I have provided an argument based on a fitness function proposed in the reference article that clearly explains the decrease. Although the main result of the reference article - that individual learning can indirectly guide evolutionary search - has been reproduced, the difference in simulation results suggests that individual learning might be eradicated over time. This decrease in usefulness of individual learning occurs in static environments where targeted genome stays fixed, as in these simulations. In more realistic scenarios, where environments are non-stationary and the targeted genome would be changing over time, individual learning would have a more stable role.


# References

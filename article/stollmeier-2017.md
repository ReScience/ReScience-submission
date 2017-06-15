---
Title: "A simple rule for the evolution of cooperation on graphs and social networks"
Author:
  - name: Frank Stollmeier
    affiliation: 1, 2
Address:
  - code:    1
    address: Network Dynamics, Max Planck Institute for Dynamics and Self-Organization (MPIDS), Am Faßberg 17, 37077 Göttingen, Germany
  - code:    2
    address: Institute for Nonlinear Dynamics, Faculty of Physics, University of Göttingen, Am Fassberg 17, 37077 Göttingen, Germany
Contact:
  - franks@nld.ds.mpg.de
Editor:
  - Name Surname
Reviewer:
  - Name Surname
  - Name Surname
Publication:
  received:  xxx,  x, xxxx
  accepted:  xxx, x, xxxx
  published: xxx, x, xxxx
  volume:    "**x**"
  issue:     "**x**"
  date:      xxx xxx
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:      
  notebook:  
Reproduction:
  - "A simple rule for the evolution of cooperation on graphs and social networks. H. Ohtsuki, C. Hauert, E. Lieberman and M.A. Nowak. Nature ***441*** (2006)."
Bibliography:
  bibliography.bib

---

# Introduction

A central question in evolutionary game theory is how cooperation can evolve in a Prisoner's Dilemma.
In an unstructured population the defectors have always a higher payoff than the cooperators, hence natural selection leads to extinction of the cooperators. One of the mechanisms that promote cooperation is network reciprocity. If the individuals do not interact with all other individuals but only with a subset of the other individuals, e.g. the adjacent individuals in a social network, then cooperators can form groups with other cooperators, which increases the benefit from other cooperators and reduces the exploitation from defectors. As a consequence, the fixation probability of cooperation increases. Cooperation is said to be favored by natural selection if this probability exceeds $1/N$, the fixation probability of a neutral mutation in a population of size $N$. 

Ohtsuki et al. [@Ohtsuki2006] showed that a simple rule indicates whether cooperation is favored. The fixation probability of cooperators exceeds $1/N$ if the condition $b/c>k$ is satisfied. The parameter $c$ is the cost an individual pays to cooperate with others, $b$ is the benefit the other individual receives from a cooperator, and $k$ is the average degree of the network. They compare this rule with fixation probabilities measured in numerical simulations. The rule is most accurate for regular networks with large populations. It is less accurate, but still valid, for networks in which the node degree varies, e.g. random or scale-free networks, and for networks with small population sizes.
Until today, the question of the critical benefit-to-cost ratio in networks attracted much interest. Further proofs of this rule, more precise rules, exact calculations and generalizations were published [@Allen2014; @Allen2017; @Chen2013; @Konno2011; @Taylor2007].

This paper presents an implementation to reproduce the numerical results of the original paper based solely on the information given in the original paper.
The computationally intensive algorithms are written in Cython. A Jupyter Notebook contains the Python scripts to produce all figures. Since accurate measurements of the fixation probability require a large number of simulations for each set of parameters, these scripts are prepared to run simulations in parallel using the module ipyparallel.


# Methods

## Overview

In order to measure the fixation probability for a certain graph type and specific parameters, we generate 1000 graphs and run 1000 simulations on each graph. 
Each simulation starts with the initialization of the strategies. All individuals are set to 'defect', except for one randomly chosen individual that is set to 'cooperate'. After the initialization the death-birth updating or imitation updating is repeated until all individuals have the same strategy. The fixation probability is the ratio of the number of simulations in which all individuals cooperate at the end and the total number of the simulations. Cooperation is expected to be favored if $b/c>k$ for the death-birth updating and if $b/c>k+2$ for imitation updating.
We use the same graph types as described in the supplementary of the original paper: circle graphs, lattice graphs, random regular graphs, random graphs and scale-free graphs.

## Death-birth update

A death-birth update consists of two steps. First, a random individual is chosen to die. Second, the adjacent individuals compete for filling the vacant site with a copy of themselves. If the adjacent cooperators have the total fitness $F_C$ and the adjacent defectors have the total fitness $F_D$, the probability that the new individual is a cooperator is $\frac{F_C}{F_C+F_D}$. The fitness of an individual is $f = 1-w+wP_X$. The parameter $w$ is the selection strength and throughout this study set to $0.01$ ("weak selection"). The payoff $P_X$ depends on the strategy $X\in \{C,D\}$ of the individual, the number of adjacent cooperators $i$ and the number of adjacent defectors $j$. A cooperator gets $P_C=bi-c(i+j)$ and a defector $P_D=bi$, where $b$ is the benefit from a cooperator and $c$ the cost of cooperation. Note that the individual that is chosen to die is included when evaluating the fitness of the adjacent individuals. In other words, the fitness is evaluated before the death step.

## Imitation updating

First, a random individual is chosen to update its strategy. Including to the total fitnesses $F_C$ and $F_D$ of the adjacent cooperators and defectors the imitation update also depends on the own fitness $f_0$ of the selected individual. With the probability $\frac{F_C}{F_C+F_D+f_0}$ it will imitate the cooperation strategy and with the probability $\frac{F_D}{F_C+F_D+f_0}$ it will imitate the defection strategy.  



# Results

In the end we are interested in how the fixation probability depends on the benefit-to-cost ratio $b/c$, but for the simulations we need to choose absolute values for the costs $c$ and the benefit $b$. Since neither of those values is given in the original paper, we first study how the choice of $c$ affects the results.
As an example, we measure the fixation probability curves for the cycle graph with $N=500$ nodes and degree $k=10$ for different choices of $c$. 
The results are shown in figure @fig:c_estimation, together with data taken from figure 2a of the original paper for comparison. 

Most importantly, the position of the intersection of the fixation probability curves with the fixation probability of a neutral mutation is almost unaffected by the choice of the absolute value $c$. Assuming that this is also true for other graph types, this means that we can choose (within a certain range) arbitrary absolute values. The curve with $c=0.125$ matches the data of the same graph type from Ohtsuki et al. We will use this value for all further simulations. 

![Fixation probability curves for a cycle graph with $N=500$ nodes and degree $k=10$. The colored dots are measured with different absolute cost values, the black dots are taken from figure 2a of the original paper.](./figures/cost_value_estimation.pdf){#fig:c_estimation}

Figure @fig:fp_db shows the fixation probability curves for the death-birth process. As expected from the theory, all curves cross the line of a neutral mutation at benefit-to-cost ratios greater than $k$. With larger networks the intersection in much closer to $b/c=k$. 
This effect is much weaker with scale-free networks than with the other network types. These observations are consistent with the original paper. The only visible difference is that the curves in original paper are much steeper.

Figure @fig:fp_im shows the fixation probability curves for imitation updating with small network sizes only, as in the original paper. All curves cross the line of a neutral mutation at benefit-to-cost ratios greater than $k+2$. Again, the difference to the original paper is that the curves in the original paper are much steeper. 

![Fixation probability curves using death-birth updating with the same graph types and parameters as in Figure 2 of the original paper. The only parameter that is different to the original paper is the absolute cost value $c$, which is $0.125$ here and unknown in the original paper. The horizontal black dashed lines indicate the fixation probability of a neutral mutation $1/N$, the vertical colored lines are at $b/c = k$.](./figures/fixation_probabilities.pdf){#fig:fp_db}

![Fixation probability curves using imitation updating with the same graph types and parameters as in Figure 4 of the supplementary of the original paper. The only parameter that is different to the original paper is the absolute cost value $c$, which is $0.125$ here and unknown in the original paper. The horizontal black dashed lines indicate the fixation probability of a neutral mutation $1/N$, the vertical colored lines are at $b/c = k+2$.](./figures/fixation_probabilities_im.pdf){#fig:fp_im}

# Conclusion

The steeper curves in the original paper are likely a result of larger absolute cost values. The test with different values showed that the absolute cost value changes the slope of the curve, but has almost no effect on the position of the intersection of the fixation probability curve with the fixation probability of a neutral mutation. For all considered graph types, network sizes and degrees the positions of the intersection are consistent with the numerical results presented in the original paper and are consistent with their theoretical results, $b/c>k$ for death-birth updating and $b/c>k+2$ for imitation updating.
In summary, the results of this implementation confirm the numerical results of the original paper. 



# References

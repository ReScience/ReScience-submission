---
Title: "How Attention Can Create Synaptic Tags for the Learning of Working Memories in Sequential Tasks"
Author:
  - name: Alexandre Frédéric
    affiliation: 1, 2
  - name: Le Masson Erwan,
    affiliation: 1, 2, 3
Address:
  - code:    1
    address: INRIA Bordeaux Sud-Ouest, Bordeaux, France
  - code:    2
    address: Institut des Maladies Neurodégénératives, Université de Bordeaux, Bordeaux, France
  - code:    3
    address: ENSEIRB-MATMECA, Talence, France
Contact:
  - frederic.alexandre@inria.fr
Editor:
  - Name Surname
Reviewer:
  - Name Surname
  - Name Surname
Publication:
  received:  Sep,  1, 2015
  accepted:  Sep, 1, 2015
  published: Sep, 1, 2015
  volume:    "**1**"
  issue:     "**1**"
  date:      Sep 2015
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:      
  notebook:  
Reproduction:
  - "How Attention Can Create Synaptic Tags for the Learning of Working Memories in Sequential Tasks,
     J. Rombouts, M. Bohte and P. Roelfsema, PLoS Computational Biology 11.3, e1004060. DOI: 10.1371/journal.pcbi.1004060"
Bibliography:
  article.bib

---

# Introduction

The reference paper [@jrombouts:2015] introduces a new reinforcement learning model using memory
units to learn sequential tasks called AuGMEnT. The results presented suggest new 
approaches in understanding the acquisition of more complex tasks, and more generic learning 
mechanisms. An implementation of the model was provided by the author which helped
verifying the correctness of intermediate computations. Python 3 was used for this replication along with NumPy
and multiprocessing libraries for speedup. Object oriented architecture is used in the module
to help factorize the code.


# Methods

## Model

The initial intention was to implement the model using an artificial neural network
simulator. The simulation tool ANNarchy [@vitay:2015] was considered for its ability to simulate rate-coded
networks. Unfortunately, the fixed order of evaluation, i.e. connexions then populations,
makes it difficult to implement back-propagation models such as AuGMEnT. It was instead
decided to write a custom script to simulate the network.\


The paper's description of the model details all the update functions and is relatively straight forward to
implement, only the initial value of $Q_{a}(t-1)$ was not provided for the first computation of $\delta$ in equation 17.
Where the author's implementation uses $Q_{a}(t)$, it was decided to use a more naive solution and set it to 0.
It might also be useful to clarify the nature of the feedback weights $w'$ in equations 14 and 16: once an action is selected,
only the feedback synapses leaving the corresponding selected Q-value unit are activated to update tags,
more precisely: $w'_{ij} = w_{ij} * z_{i}$. The model can have specific feedback synapses as well but the simpler method is to
use the feed-forward synapses' weights.\


To offer some discussion about the model and its limits, the first point to bring forward would be its artificial time management.
The extreme discretization of time and explicit signals such as trial begin and end make it difficult to consider real-time simulation
or even realistic environments implementations. Another point worth noting is the possible ambiguity of memory traces. Because the traces in memory units are the sums of
changes in input, there exist some sequences of inputs the model would be incapable of distinguishing. For example, the
sequences `((0, 0), (1, 0), (1, 1))` and `((0, 0), (0, 1), (1, 1))` have the same memory traces `(1, 1)`.

## Tasks

The descriptions of the tasks used to test the network are somewhat minimal and it is sometimes
necessary to read the refer to the original experiments' papers [@gottlieb:1999] and [@yang:2007] for clarifications.
In this section, some details of implementation will be exposed.\


For the fixation tasks , i.e. saccade/anti-saccade and probabilistic decision making, once the network has
fixated the point, it has to maintain fixating until the "go" signal where the screen turns off, otherwise the experiment is failed.
Moreover, during the "go" phase, the gaze can only be chosen once, if it is not the target, the trial fails.
At the end of the trial, a blank screen is shown before exiting.\


The shaping strategy for the probabilistic decision making task consists in increasing gradually the difficulty
of the task. The table 3 in the article describes all 8 levels of difficulty. The column "# Input Symbols" is
the size of the subset of shapes. The network is not confronted to all shapes immediately:
first, the two infinite weight shapes are used, then shapes with the smallest absolute weights are added as the difficulty increases.
The column "Sequence Length" is the number of shapes shown during a trial. The more shapes there are on screen, the more difficult
it is to determine which target should be chosen.
A number of settings are randomized, such as which shapes should appear, their order of apparition, but also their locations around
the fixation point. The first shape can appear in any of the 4 locations, the second in any of the remaining 3 locations, etc.
Finally, the triangle and heptagon, shapes with infinite weights, cancel each other in the computation of the total weight.


# Results

Only the saccade/anti-saccade task and the probabilistic decision making task were implemented.
The results observed, presented in table @tbl:results, are comparable to the original article's.
The differences could come from undocumented changes in the experiments' protocols. Qualitative results,
for example the use of shaping strategy to obtain better performances, are confirmed by this replication.


Task                Success in [@jrombouts:2015] Success Convergence in [@jrombouts:2015] Convergence
------------------- ---------------------------- ------- -------------------------------- --------------
Saccade w/ shaping  99.45%                       85.0%   4100 trials                      3800 trials
Saccade w/o shaping 76.41%                       64.0%   N/A                              5433 trials
Probabilistic       99.0%                        100%    55234 trials                     70762.5 trials
------------------- ---------------------------- ------- -------------------------------- --------------
Table: Results {#tbl:results}


# Conclusion

The results obtained are comparable to those announced in the article. Ambiguities in the
experiments' descriptions could be the cause for some variations, but do not contradict the
article's conclusion.


# References

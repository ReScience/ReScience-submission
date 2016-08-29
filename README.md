### Reference Implementation of [1]

The source code is written in Python (version 3.5.1) using Numpy (version
1.10.4), Scipy (version 0.17.0) and Matplotlib (version 1.5.1). The source 
code is distributed under the GPLv3 and a copy of the License is provided 
along with the source code. The repository is organized as follows:
```
├── article
│   ├── article_bak.tex
│   ├── article.pdf
│   ├── article.tex
│   ├── crossref.yaml
│   ├── detorakis-2016.bcf
│   ├── detorakis-2016.bib
│   ├── detorakis-2016.md
│   ├── detorakis-2016.pdf
│   ├── detorakis-2016.run.xml
│   ├── detorakis-2016.tex
│   ├── figs
│   │   ├── const.pdf
│   │   ├── Figure1.pdf
│   │   ├── Figure2.pdf
│   │   ├── Figure3.pdf
│   │   ├── Figure4.pdf
│   │   ├── Figure5.pdf
│   │   └── pulse.pdf
│   ├── LICENSE
│   ├── pandoc-template.tex
│   ├── README.md
│   ├── rescience-logo.pdf
│   └── rescience-template.tex
├── code
│   ├── LICENSE
│   ├── neuron_model.py
│   ├── params
│   │   ├── params_figure1.cfg
│   │   ├── params_figure2.cfg
│   │   ├── params_figure3.cfg
│   │   ├── params_figure4b.cfg
│   │   ├── params_figure4.cfg
│   │   ├── params_figure5a.cfg
│   │   ├── params_figure5b.cfg
│   │   └── params_figure5c.cfg
│   ├── plot_figures.py
│   ├── README.md
│   ├── run_all.py
│   ├── simulations.py
│   ├── stimulus.npy
│   └── tools.py
├── data
│   ├── data1.dat
│   ├── data.dat
│   ├── LICENSE
│   └── README.md
├── lalal
├── LICENSE
├── notebook
│   └── README.md
└── README.md
```

* **article** -> This folder contains the accompanying text.
* **code** -> Here is the source code along with all the plotting functions.
* **data** -> This folder is normally empty, but it is used by the source 
  code for storing the results of the simulations. 
* **notebook** -> There is no available notebook for this implementation. 



[1] : *Multiple dynamical modes of thalamic relay neurons: rhythmic
bursting and intermittent phase-locking*, Wang, X-J, Neuroscience,
59(1), pg.21-31, 1994.

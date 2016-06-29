## Robust timing and motor patterns by taming chaos in recurrent neural networks

**Author:** Julien Vitay <julien.vitay@informatik.tu-chemnitz.de>

*Professorship for Artificial Intelligence, Department of Computer Science, Chemnitz University of Technology, D-09107 Chemnitz, Germany.*

A reference implementation of:

    Laje, R. and Buonomano, D.V. (2013). Robust timing and motor patterns by taming chaos in recurrent neural networks. Nat Neurosci. 2013 Jul;16(7):925-33. <doi:10.1038/nn.3405>

The original article and the associated data/code can be found online on Pubmed: <http://www.ncbi.nlm.nih.gov/pubmed/23708144>

### Dependencies

The standard scientific Python stack is required:

* Python 2.7 or >= 3.4
* Numpy 1.10 (lower versions may work but not tested)
* Scipy 0.17
* Matplotlib 1.3

### Data

The handwriting patterns for Fig. 2 are available on PMC (<http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3753043>), but the copyright holder is the Nature Publishing Group and no free license is provided, so it cannot be included in this repository. In order to reproduce Fig. 2, one has to download the provided data to obtain a `.mat` file.

The data is located at <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3753043/bin/NIHMS472497-supplement-3.zip>. This zip file should then be decompressed in the `data/` older, so that the file `DAC_handwriting_output_targets.mat` lies there.

In `data/` is provided a `get_handwriting.sh` script for Linux/Mac OS users that automatically performs these steps:

```{.bash}
wget http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3753043/bin/NIHMS472497-supplement-3.zip -O data.zip
unzip data.zip DAC_handwriting_output_targets.mat
```

### Model

The model is implemented by a class `RecurrentNetwork` in the file `code/RecurrentNetwork.py`. The scripts `code/Fig1.py`, `code/Fig2.py` and `code/Fig3.py` allow to reproduce the corresponding figures of the manuscript.

As the script for Figure 3 takes 3 days of computation on a standard computer, we provide the script `code/Fig3-Load.py` that only produces the figure, based on recoreded data stored in `data/timingcapacity.npz`.

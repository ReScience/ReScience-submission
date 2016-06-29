#!/bin/bash
# Script allowing to download the data needed for Fig. 2 of:
#
#   Laje, R. and Buonomano, D.V. (2013). Robust timing and motor patterns by taming chaos in recurrent neural networks. Nat Neurosci.
#
# Author: Julien Vitay (julien.vitay@informatik.tu-chemnitz.de)
# Licence: MIT

# Download the archive from PMC
wget http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3753043/bin/NIHMS472497-supplement-3.zip -O data.zip

# Unzip only the MAT file
unzip data.zip DAC_handwriting_output_targets.mat

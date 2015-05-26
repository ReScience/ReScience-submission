# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2015, Meropi Topalidou
# Distributed under the (new) BSD License.
#
# Contributors: Meropi Topalidou (Meropi.Topalidou@inria.fr)
#				Nicolas Rougier  (Nicolas.Rougier@inria.fr)
# -----------------------------------------------------------------------------
# References:
#
# * Interaction between cognitive and motor cortico-basal ganglia loops during
#   decision making: a computational study. M. Guthrie, A. Leblois, A. Garenne,
#   and T. Boraud. Journal of Neurophysiology, 109:3025–3040, 2013.
# -----------------------------------------------------------------------------

# Population size
n = 4
n_trials = 120
n_learning_trials = 120
n_testing_trials  = 120
simulations = 250
# --- Time ---
ms           = 0.001
duration	 = int(3.0/ms)
dt           = 1*ms
tau          = 10*ms

# --- Learning ---
alpha_CUE  = 0.050
alpha_LTP  = 0.01
alpha_LTD  = 0.005

# --- Sigmoid ---
Vmin = 0
Vmax = 20
Vh   = 16
Vc   = 3

# --- Model ---
decision_threshold = 40
noise              = 0.001
CTX_rest   =  -3.0
STR_rest   =   0.0
STN_rest   = -10.0
GPI_rest   =  10.0
THL_rest   = -40.0

# Noise level (%)
Cortex_N   =   0.01
Striatum_N =   0.001
STN_N      =   0.001
GPi_N      =   0.03
Thalamus_N =   0.001

# --- Cues & Rewards ---
V_cue   = 7
rewards_Guthrie = 3/3.,2/3.,1/3.,0/3.
rewards_Guthrie_reverse_all = 0/3.,1/3.,2/3.,3/3.
rewards_Guthrie_reverse_middle = 3/3.,1/3.,2/3.,0/3
rewards_Piron  = 0.75, 0.25, 0.75, 0.25
rewards_Piron_reverse  = 0.25, 0.75, 0.25, 0.75

# -- Weight ---
Wmin  = 0.25
Wmax  = 0.75
gains = { "CTX.cog -> STR.cog" : +1.0,#0.0,#
          "CTX.mot -> STR.mot" : +1.0,
          "CTX.ass -> STR.ass" : +1.0,
          "CTX.cog -> STR.ass" : +0.2,
          "CTX.mot -> STR.ass" : +0.2,
          "CTX.cog -> STN.cog" : +1.0,
          "CTX.mot -> STN.mot" : +1.0,
          "STR.cog -> GPI.cog" : -2.0,
          "STR.mot -> GPI.mot" : -2.0,
          "STR.ass -> GPI.cog" : -2.0,
          "STR.ass -> GPI.mot" : -2.0,
          "STN.cog -> GPI.cog" : +1.0,
          "STN.mot -> GPI.mot" : +1.0,
          "GPI.cog -> THL.cog" : -0.5,
          "GPI.mot -> THL.mot" : -0.5,

          "THL.cog -> CTX.cog" : +1.0,
          "THL.mot -> CTX.mot" : +1.0,
          "CTX.cog -> THL.cog" : +0.4,
          "CTX.mot -> THL.mot" : +0.4,

 }
dtype = [ ("CTX", [("mot", float, 4), ("cog", float, 4), ("ass", float, 16)]),
          ("STR", [("mot", float, 4), ("cog", float, 4), ("ass", float, 16)]),
          ("GPI", [("mot", float, 4), ("cog", float, 4)]),
          ("THL", [("mot", float, 4), ("cog", float, 4)]),
          ("STN", [("mot", float, 4), ("cog", float, 4)])]

threshold  = 40

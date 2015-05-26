# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2015, Meropi Topalidou
# Distributed under the (new) BSD License.
#
# Contributors: Meropi Topalidou (Meropi.Topalidou@inria.fr)
#				Nicolas Rougier  (Nicolas.Rougier@inria.fr)
# -----------------------------------------------------------------------------

# Evolution of single trial with Guthrie protocol
# -----------------------------------------------------------------------------

if __name__ == "__main__":
	from model import *
	from display import *
	from trial import *

	cues_pres = input('\nDo you want to present cues?\nChoose 1 for True or 0 for False\n')
	reset(protocol = 'Guthrie')
	# Make GPi lesion
	# connections["GPI.cog -> THL.cog"].active = False
	# connections["GPI.mot -> THL.mot"].active = False
	histor, time = trial(hist = True, debugging = True, protocol = 'Guthrie', cues_pres = cues_pres, wholeFig = True)
	if 1: display_ctx(histor, 3.0)
	if 0: display_ctx(histor, 3.0, "single-trial-NoBG.pdf")
	if 1: display_all(histor, 3.0)#, "single-trial-all.pdf")

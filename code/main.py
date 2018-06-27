# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# This script creates all figures. if run without the "--new True" flag, 
# Figures 5 and 6 will be created using the saved data in the data folder. 
# Usage: python main.py 
# Usage: python main.py --new True
# To change the number of runs for figures 5 and 6, change the reps_per_combination
# for figure 5 and reps_per_resolution for figure 6 to a lower number.
# -----------------------------------------------------------------------------

from figure_1 import fig_1
from figure_2 import fig_2
from figure_5 import fig_5_saved, fig_5_new
from figure_6 import fig_6_saved, fig_6_new
from run_simulation import run_sim
from figure_7CD import fig_7CD
from figure_8 import fig_8
from write_param_inputfile import write_para_file
import matplotlib.pyplot as plt
import argparse
import time

# ### get command line options
ap = argparse.ArgumentParser()
ap.add_argument("-n", "--new", required=False, default='False',
                help="whether to run calculations for figures 5 and 6 (True) or use saved values (False, default)")
args = vars(ap.parse_args())

# ### plot figures from paper
fig1 = fig_1()

fig2 = fig_2()

print('%s Preparing Figures 3, 4, 7AB' % time.strftime('%H:%M'))
fig3, fig4, fig7AB, fig9, fig10 = run_sim(1)

if args['new'] == 'False':
    fig5 = fig_5_saved()
    fig6 = fig_6_saved()

elif args['new'] == 'True':
    write_para_file()

    print('%s Preparing Figure 5' % time.strftime('%H:%M'))
    reps_per_combination = 10
    fig5 = fig_5_new(reps_per_combination)

    print('%s Preparing Figure 6' % time.strftime('%H:%M'))
    reps_per_resolution = 10
    fig6 = fig_6_new(reps_per_resolution)

else:
    'input for --new was not recognised. Please try again'

fig7CD = fig_7CD(28)

fig8 = fig_8()

plt.show()

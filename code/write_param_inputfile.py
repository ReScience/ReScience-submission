# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Writes the parameter file that stores the combinations of parameters used for 
# Figure 5.
# This script can be run independently.
# -----------------------------------------------------------------------------

from brian2 import *

'''
Parameter combinations as saved in param.npy

row     win     jit     # pat neur      pat freq        % deleted
0       0.475   1       1000            0.25            0.0

1       0.275   1       1000            0.25            0.0
2       0.325   1       1000            0.25            0.0
3       0.375   1       1000            0.25            0.0
4       0.425   1       1000            0.25            0.0

5       0.475   0       1000            0.25            0.0
6       0.475   2       1000            0.25            0.0
7       0.475   3       1000            0.25            0.0
8       0.475   4       1000            0.25            0.0
9       0.475   5       1000            0.25            0.0
10      0.475   6       1000            0.25            0.0

11      0.190   1       400             0.25            0.0
12      0.285   1       600             0.25            0.0
13      0.380   1       800             0.25            0.0
14      0.570   1       1200            0.25            0.0

15      0.475   1       1000            0.05            0.0
16      0.475   1       1000            0.10            0.0
17      0.475   1       1000            0.15            0.0
18      0.475   1       1000            0.50            0.0

19      0.4275  1       1000            0.25            0.1
20      0.380   1       1000            0.25            0.2
21      0.3325  1       1000            0.25            0.3

'''


def write_para_file():
    r_ = array([4])
    w_ = array([0.275, 0.325, 0.375, 0.425])
    j_ = array([0, 2, 3, 4, 5, 6])
    n_ = array([400, 600, 800, 1200])
    f_ = array([0.05, 0.10, 0.15, 0.50])
    d_ = array([0.10, 0.20, 0.30])
    combinations = len(r_) + len(w_) + len(j_) + len(n_) + len(f_) + len(d_)
    stand = array([4, 0.475, 1, 1000, 0.25, 0, 500])

    para = empty((combinations, 7))
    for i in range(combinations):
        para[i] = stand
    count = 0
    for x in r_:
        para[count][0] = x
        count += 1
    for x in w_:
        para[count][1] = x
        count += 1
    for x in j_:
        para[count][2] = x
        count += 1
    for x in n_:
        para[count][3] = x
        para[count][6] = 0.5 * (1 - stand[5]) * x
        para[count][1] = 1.9 * 0.5 * (1 - stand[5]) * x / 2000
        count += 1
    for x in f_:
        para[count][4] = x
        count += 1
    for x in d_:
        para[count][5] = x
        para[count][6] = 0.5 * (1 - x) * stand[3]
        para[count][1] = 1.9 * 0.5 * (1 - x) * stand[3] / 2000
        count += 1
    # print(para)
    save('../data/para.npy', para)

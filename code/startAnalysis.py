"""
Python script to run all experiments after each other to reproduce the presented data
in Clopath et al. (2010)
"""

import os

def start():
    print('Start to reproduce the data in Clopath et al. (2010)')
    os.system('python Fig1_clamp.py --clean')
    os.system('python Fig1_window.py --clean')
    os.system('python Fig1_pairing.py --clean')
    os.system('python Fig2_burst.py --clean')
    os.system('python Fig3_rateCode.py --clean')
    os.system('python Fig3_temporalCode.py --clean')
    os.system('python Fig3_rateCode_stand.py --clean')
    os.system('python Fig3_temporalCode_stand.py --clean')
    os.system('python Fig4_stableW.py --clean')
    if os.path.isfile('IMAGES.mat'):
        os.system('python Fig4_RF.py --clean')
    else:
        print('No IMAGES.mat found, please download the file from: https://www.rctn.org/bruno/sparsenet/IMAGES.mat')
    print('Done with running all the scripts.')
if __name__ == "__main__":
    start()

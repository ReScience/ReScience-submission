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
    os.system('python Fig4_stableW.py --clean')
    os.system('python Fig4_RF.py --clean')

    print('Done with running all the scripts.')
if __name__ == "__main__":
    start()

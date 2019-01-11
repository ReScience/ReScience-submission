#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Cressot Loic
    ISIR - CNRS / Sorbonne Université
    12/2018
""" 

"""
    gen_test : simple test to confirm the good working of data generation from gym environment gym_round_bot
    Note : This script should not be called manually but rather should be called by the bash script gen_test.sh
"""

import os, sys
# add parent folder to path
sys.path.append(sys.path[0] + '/..')
import tools
import argparse

parser = argparse.ArgumentParser(description=None)
parser.add_argument('-d','--datafile', default='', help='Select data file to load', required=True)
parser.add_argument('-rt','--recordto', default='', help='Select folder to record data', required=True)
parser.add_argument('-n','--name', default='', help='name for the image', required=True)
args = parser.parse_args()

data = tools.load_data(args.datafile)  
tools.plot_observations(  data['observations'], 
                          name="Test on generation of observations",
                          plot_name=args.name,
                          recordto=args.recordto,
                          display=False,
                        )

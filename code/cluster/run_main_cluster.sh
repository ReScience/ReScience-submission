#!/bin/bash

#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N rescience

module load anaconda3/personal

python /THIS/PATH/NEEDS/TO/BE/CHANGED/main_cluster.py runid=$PBS_ARRAY_INDEX

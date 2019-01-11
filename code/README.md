## Introduction

***[Re] Learning state representations with robotic priors***  
Loic Cressot, Alexandre Coninx, Astrid Merckling, Nicolas Perrin, and Stéphane Doncieux  
*Sorbonne Université, CNRS, Institut des Systèmes Intelligents et de Robotique, ISIR, F-75005 Paris, France*  
  
**A reference implementation of:**  
*Learning state representations with robotic priors, Jonschkowski, R. and Brock, O., 2015. Autonomous Robots, 39(3), pp.407-428*


## Requirements

+ **python version:** python3

+ **system:** Linux or MacOS

+ **required packages along with a working version** 
(all these packages can be installed with pip):
	+ *gym:* 0.10.5
	+ *numpy:* 1.14.2
	+ *pillow:* 5.3.0
	+ *tensorflow:* 1.5.0
	+ *keras:* 2.1.5
	+ *tqdm:* 4.23.4
	+ *scikit-learn:* 0.19.1
	+ *matplotlib:* 2.1.2

+ **additional package:**
gym_round_bot 0.0.1 (requires: scipy 1.0.1,  pyglet 1.3.1)  
This is a gym environment we coded. It is needed to generate training sets and to evaluate policies.
It can be downloaded at https://github.com/Lcressot/gym-round_bot and should be easily installable on UNIX via terminal with ```pip install -e .``` once inside the directory.

+ **installation advice on Ubuntu:**  
Some libraries installation issues may be encountered on Linux. We advice the creation of a tensorflow python 3 virtual environment with Anaconda, then to install all possible libraires with the conda installer (except for the gym librairies which have to be installed with pip3).  
Furthermore, for some Linux configurations (e.g. Ubuntu 18.04 and NVIDIA driver 390.77), it seems that the display and pixel value acquisition are not working properly. We could not find the problem which has probably something to do with the graphic driver, and suggest to run beforehand the test script gym_env/gen_test.sh and look at the generated image gym_env/temp_test/observations.png which should be coloured and not dark or grey.

## Code organization

This repository is organized as follows :
+ main.py is the used to reproduce the experiments of the article.
+ jonschkowski_priors.py implements the state representation learning method presented in the article.
+ fqiteration.py implements the fitted q iteration method presented in the article.
+ tools.py contains usefull functions for loading data, and for computing and plotting representations.
+ gym_env is a folder containing some code for generating new traning data with our gym_round_bot OpenAI gym environment. The script for generating this training data is main_agent_vs_env.py.

## Usage

**training set generation:**
Use the gym_round_bot environment and the script gym_env/main_agent_vs_env.py to generate a training set:
```Bash
python3 gym_env/main_agent_vs_env.py --help
```
**Note:** you can check the good functionning of the data generation by running the script gym_env/gen_test.sh and verifying the generated image gym_env/temp_test/observations.png

More info on https://github.com/Lcressot/gym-round_bot

**state representation learning:**
Learn a state representation on a training set:
```Bash
python3 main.py --help
```

**fitted Q iteration on learned representation:**
Evaluate a state representation at different training steps with fitted Q iteration:
```Bash
python3 main.py -ql --help
```

## Experiment reproduction

Here we present the command line parameters we used to produce the experiment of the article, both for the training set generation and the state representation learning.

### Simple navigation task 

training set generation:
```Bash
python3 gym_env/main_agent_vs_env.py -c XZ --world_name square --texture colours -xzt 2 2 1 -ms 50 -f 150 -nr 0.1 -obs 16 16 -ml -75 75 -nep 100 -rt ../data/ --world_size 45 45 --speed 3 --name simple
```

state representation learning:
```Bash
python3 main.py -r ./results/ -ne 10 -lr 0.001 -reg 0.03 -sd 2 -bs 256 -dis -trd ../data/simple_XZ_2x2x1_16x16_50x100_s3_f150_square_ws45x45_colours_noise0.10_ml2.npz 
```

fitted Q iteration on learned representation:
```Bash
python3 main.py -r ./results/ -ne 25 -lr 0.0001 -reg 0.03 -sd 2 -bs 256 -trd ../data/simple_XZ_2x2x1_16x16_50x100_s3_f150_square_ws45x45_colours_noise0.10_ml2.npz -ql
```


### Invariance to perspective

training set generation:
```Bash
python3 gym_env/main_agent_vs_env.py -c XZ --world_name square --texture colours -xzt 2 2 1 -ms 50 -f 150 -nr 0.1 -obs 16 16 -nep 100 -rt ../data/ --world_size 45 45 --speed 3 -agp -rd 6 --name topdown
```

state representation learning:
```Bash
python3 main.py -r ./results/ -ne 20 -lr 0.001 -reg 0.0001 -sd 2 -bs 256 -dis -trd ../data/topdown_XZ_2x2x1_16x16_50x100_s3_f65_square_ws45x45_gp-True_colours_noise0.10.npz 
```


### Ignoring distractors

training set generation:
```Bash
python3 gym_env/main_agent_vs_env.py -c XZ --world_name square --texture colours -xzt 2 2 1 -ms 50 -f 150 -nr 0.1 -obs 16 16 -ml -75 75 -nep 100 -rt ../data/ --world_size 45 45 --speed 3 --distractors --name distractors
```

state representation learning:
```Bash
python3 main.py -r ./results/ -ne 70 -lr 0.0001 -reg 0.1 -sd 2 -bs 256 -dis -trd ../data/distractors_XZ_2x2x1_16x16_50x100_s3_f150_square_ws45x45_colours_noise0.10_ml2_dis.npz
```

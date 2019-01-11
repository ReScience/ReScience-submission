# Data repository

**Note :** you can first check the good functionning of the data generation by running the script gym_env/gen_test.sh and verifying the generated image gym_env/temp_test/observations.png

## Simple navigation task

**data file name :** simple_XZ_2x2x1_16x16_50x100_s3_f150_square_ws45x45_colours_noise0.10_ml2.npz

**data file generation command :**
```Bash
python3 gym_env/main_agent_vs_env.py -c XZ --world_name square --texture colours -xzt 2 2 1 -ms 50 -f 150 -nr 0.1 -obs 16 16 -ml -75 75 -nep 100 -rt ../data/ --world_size 45 45 --speed 3 --name simple
```

**data parameters :**
+ controller (-c) : XZ. The robot moves on the OXZ plane always facing in the same direction.
+ world name (--world_name) : square. Choose the square world for this experiment.
+ texture (--texture) : colours. Apply uniform colours texture to the blocks of this world.
+ xzt range (-xzt) : 2x2x1. Range for X and Z dispclacements and Theta rotation (unused for XZ controller). With a speed of 3, the displacements can be of [-6, -3, 0, 3, 6] units (2 x direction_range + 1 possibilities ) for each direction X and Z.
+ observation size (-obs) : 16 16. Observations are 16x16 square images with RGB colours
+ max step (-ms) : 50. Number of steps of each episode.
+ focal (-f) : 150. Focal length of frames. Here we also use multi-view rendering with two frames being squeezed and connected together horizontally.
+ noise ratio (-nr) : 0.1. The ratio to compute the standard variation of displacement additive gaussian noise from displacement's amplitude.
+ multi-view (-ml) : -75 75. As 75 is equal to 150/2 and the focal being 150°, we set the two angles to theses values for rendering relatively to the robot and obtain a final 300° image.
+ number of episodes (-nep) : 100.
+ record to (-rt) : ../data/
+ world size (--world_size) : 45 45. Width and depth of the square world in units, as in the article.
+ speed 3 (--speed) : speed in units, along one direction, of the robot. In combination with the two displacement directions and the speed range factor, the speed on each direction can be within [-6, -3, 0, 3, 6].
+ name (--name) : simple. Just a suffix name to be added in front of all parameters in the file name.



## Invariance to perspective

**data file name :** topdown_XZ_2x2x1_16x16_50x100_s3_f65_square_ws45x45_gp-True_colours_noise0.10.npz

**data file generation command :**  
```Bash
python3 gym_env/main_agent_vs_env.py -c XZ --world_name square --texture colours -xzt 2 2 1 -ms 50 -f 150 -nr 0.1 -obs 16 16 -nep 100 -rt ../data/ --world_size 45 45 --speed 3 --name new -agp -rd 6 --name topdown
```

**data parameters (different from those of the simple navigation task) :**
+ name (--name) : topdown. Just a suffix name to be added in front of all parameters in the file name.
+ multi-view (-ml) : not needed here because the rendering is made from a global point of view (top down)
+ automatic global point of view (-agp) : computes automatically the global point of view to be able to see the whole setup at once from top down camera.
+ radius of the robot (-rd) : 6. Set a bigger radius (3 times bigger) to the robot than in the simple navigation task to clearly see it with top down camera.



## Ignoring distractors

**data file name :** distractors_XZ_2x2x1_16x16_50x100_s3_f150_square_ws45x45_colours_noise0.10_ml2_dis.npz

**data file generation command :**  
```Bash
python3 gym_env/main_agent_vs_env.py -c XZ --world_name square --texture colours -xzt 2 2 1 -ms 50 -f 150 -nr 0.1 -obs 16 16 -ml -75 75 -nep 100 -rt ../data/ --world_size 45 45 --speed 3 --name dis --distractors --name distractors
```
**data parameters (different from those of the simple navigation task) :**
+ name (--name) : distractors. Just a suffix name to be added in front of all parameters in the file name.
+ add visual distractors (-dis) : visual distractors, i.e randomly moving coloured rectangles, are added. They move randomly on the surface of the walls and flour.

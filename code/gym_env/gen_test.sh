# gen_test : simple test script to confirm the good working of data generation from gym environment gym_round_bot

# get working directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# create or clean temp directory :
mkdir $DIR/temp_test
rm $DIR/temp_test/*

# generate data :
python3 $DIR/main_agent_vs_env.py -c XZ --world_name square --texture colours -xzt 2 2 1 -ms 50 -f 150 -nr 0.1 -obs 16 16 -ml -75 75 -nep 10 --world_size 45 45 --speed 3 --name test1 -rt $DIR/temp_test/
python3 $DIR/main_agent_vs_env.py -c XZ --world_name square --texture colours -xzt 2 2 1 -ms 50 -f 150 -nr 0.1 -obs 16 16 -nep 10 --world_size 45 45 --speed 3 --name test2 -rt $DIR/temp_test/

# display it on screen :
#chmod a+rw `$DIR/temp_test/test*.npz`
#python3 $DIR/gen_test.py --datafile `$DIR/temp_test/test1*.npz` --recordto $DIR/temp_test/ --name test1
python3 $DIR/gen_test.py --datafile $DIR/temp_test/test1_XZ_2x2x1_16x16_50x10_s3_f150_square_ws45x45_colours_noise0.10_ml2.npz --recordto $DIR/temp_test/ --name test1
#python3 $DIR/gen_test.py --datafile `$DIR/temp_test/test2*.npz` --recordto $DIR/temp_test/ --name test2
python3 $DIR/gen_test.py --datafile $DIR/temp_test/test2_XZ_2x2x1_16x16_50x10_s3_f150_square_ws45x45_colours_noise0.10.npz --recordto $DIR/temp_test/ --name test2

# clear everything :
echo ''
echo ''
echo You can now verify that generation worf by looking at image : image in $DIR/temp_test/observations.png
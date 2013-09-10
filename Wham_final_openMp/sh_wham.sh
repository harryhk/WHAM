#!/bin/bash


## generate pos_F.xvg 
#
#id=$1
#
#sed '1 d ' ../cpp_points/free$id.xvg  | awk '{print $2}' > t ; paste pos.xvg t > pos_F.xvg
#rm t 
#




# generate data in ./data folder. 

data_dir="../../bias1"

trj_dir="../../dist_arg_com1/match_bias/"


win_n=37

# ith of last j blocks

tail_n=$1 

block=$2  # 20ns block . seems too many, let's pick every 2 points. 

tail_b=$(( $tail_n * $block )) 

for(( i=0; i<=$win_n; i++  )); do 
	tail -n $tail_b $data_dir/bias_$i.xvg | head -n $block  > ./data/win_$i.xvg
	tail -n $tail_b $trj_dir/win_match_$i.xvg | head -n $block | awk ' {print $5}' > ./data/traj_$i.xvg 

done 



./wham.exe pullx.xvg -0.4 2.0 240 0.001 1000000 323 pmf_z${block}_$tail_n.xvg pos.xvg


#./wham.exe pullx.xvg 0 3.6 360 0.001 100000 323 free.xvg

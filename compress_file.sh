#! /usr/bin/bash

smpDIR=$1

ls $smpDIR/*.xls | while read file
do
	gzip $file
done

# ls $smpDIR/X1_1.0_Y1_1.0/MR_G1_X_Y/*.xls | while read file
# do
# 	gzip $file
# done


# MixRatioList="0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
# for ratio in $MixRatioList
# do
# 	ls $smpDIR/X1_${ratio}_Y1_$ratio/*.xls | while read file
# 	do
# 		gzip $file
# 	done
# done

##  Before:


## After:



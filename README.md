# HBBCHE002_FinalYear_Project_TrafficRadar
This repository contains all the MATLAB code produced for a high-resolution radar detection processor for multi-vehicle traffic scenarios

# Requirements:
- All .m and .mlx files were coded in MATLAB R2022a
- This code was created specifically for high-resolution datasets provided by the CSIR. These files would be required to run these scripts
- Most files are provided in both a .m or .mlx format 

# Steps to follow:
1. The first step is to run all the ground truth scripts. These generate .mat ground truths that are required in the other scripts. It should be noted that all these scripts generate over 40 plots each. Additionaly, its required that "groundTruth_creation_DS6p1", "groundTruth_creation_DS6p2" and "groundTruth_creation_DS6p3" are run before "groundTruth_creation_DS6"
2. After this, any other script can be run. It is recommended to run the CFAR algorithm processor scripts first.
3. It should be noted that all scripts, when run, require you to choose the dtaset you wish to run the processor on. It is thus ideal to save the dataset .mat files in the same folder

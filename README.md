# HBBCHE002_FinalYear_Project_TrafficRadar
This repository contains all the MATLAB code produced for a high-resolution radar detection processor for multi-vehicle traffic scenarios

# Requirements:
- All .m and .mlx files were coded in MATLAB R2022a
- This code was created specifically for high-resolution datasets provided by the CSIR. These files would be required to run these scripts
- Most files are provided in both a .m or .mlx format 

# Steps to follow:
1. The first step is to run all the ground truth scripts. These generate .mat ground truths that are required in the other scripts. It should be noted that all these scripts generate over 40 plots each. Additionaly, its required that "groundTruth_creation_DS6p1", "groundTruth_creation_DS6p2" and "groundTruth_creation_DS6p3" are run before "groundTruth_creation_DS6"
2. After this, any other script can be run. It is recommended to run the CFAR algorithm processor scripts first.
3. It should be noted that all scripts, when run, require you to choose the CSIR dataset you wish to run the processor on. It is thus ideal to save the dataset .mat files in the same folder

# Scripts provided:
- ground truths: scripts that generate ground truth target locations for each of the six datasets
- CA-, OS- and TM-CFAR processor scripts: these run each of the mentioned CFAR algorithms on a chosen dataset. These scripts include time domain processing, frequency domain process and an AND gate combination of the two domain processing. It contains a verification step of testing the algorithms against simulated datasets. Additionally, it runs a detection performance test on the determined detection locations by each algorithm, using the ground truth .mat files. This produces the number of false alarms, detection accuracy and number of missed detections of the algorithm
- SNR estimate script: this script calculates the SNR of each script and plots the required SNR to achieve levels of detection probability
- Threshold plotter: this script plots the calculated threshold of each algorithm alongside the dataset for a specific range line of a chosen dataset. The range line can be changed as required.  

README
11/2015 For MENG project on Matlab GPU - Testing setL_nord

Please run the following commands after initiate Matlab:
matlabpool local 12
addpath /home/NetID/yourdirectory/hetero00.71-yg326.asym.1st_do_V
addpath /home/NetID/yourdirectory/FredSigworth/EMIODist
addpath /home/NetID/yourdirectory/Asym.test.for.GPU.group
AsymObj_lmax2pmax5_create_example.m
inst_Neta2_AsymObj_rule848_Nv500_fw()

% function setL_nord will be called when you run inst_Neta2_AsymObj_rule848_Nv500_fw()
# dicomsegmentation
Python scripts for automated segmentation of DICOM files. These scripts take DICOM files, fuse them, and then segment tissue. It also calculates the standard 10-20 electrode positions for skull data, based on image fusion. 

Python:
bloodHunt.py: This file takes DICOM files and seeks out blood vessels.

MATLAB: 
testLoad.m: The master script, currently hard coded for a particular file name. Able to load files, fit 10-20 electrodes to the fused solid, and then plot the results of stimulation power in 3D. 

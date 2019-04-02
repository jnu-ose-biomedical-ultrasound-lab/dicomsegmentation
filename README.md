# dicomsegmentation
Python scripts for automated segmentation of DICOM files. These scripts take DICOM files, fuse them, and then segment tissue. It also calculates the standard 10-20 electrode positions for skull data, based on image fusion. 

Python:
bloodHunt.py: This file takes DICOM files and seeks out blood vessels.

MATLAB: 
calculateParameters.m: Calculate basic spatial parameters for a fused, 3D matrix of imported DICOM files.
positionElectrodes.m: Position 10-20 electrodes based on orientation of the head. 

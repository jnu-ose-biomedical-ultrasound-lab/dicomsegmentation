%% remove remnants of earlier iterations
close all;
clear;
clc;

%% load fused dicom model
load('rawdicom.mat')
dicom=array;

%% set initial parameters
electrodePlace=calculateParameters(dicom);
electrodePlace=positionElectrodes(electrodePlace);
cz=electrodePlace.cz;

%% electrode angle calculation
electrodePlace.angles.c3=angleBetween(cz,electrodePlace.c3);
electrodePlace.angles.c4=angleBetween(cz,electrodePlace.c4);

electrodePlace.angles.t3=angleBetween(cz,electrodePlace.t3);
electrodePlace.angles.t4=angleBetween(cz,electrodePlace.t4);

electrodePlace.angles.oz=angleBetween(cz,electrodePlace.oz);
electrodePlace.angles.pz=angleBetween(cz,electrodePlace.pz);

electrodePlace.angles.fz=angleBetween(cz,electrodePlace.fz);
electrodePlace.angles.fpz=angleBetween(cz,electrodePlace.fpz);

electrodePlace.angles.f3=angleBetween(cz,electrodePlace.f3);
electrodePlace.angles.f4=angleBetween(cz,electrodePlace.f4);

electrodePlace.angles.p3=angleBetween(cz,electrodePlace.p3);
electrodePlace.angles.p4=angleBetween(cz,electrodePlace.p4);


%% calculation of vectors 
[electrodePlace.vectors.c3,electrodePlace.startPoints.c3]=drillSpin(electrodePlace.angles.c3,cz,electrodePlace.c3,dicom);
[electrodePlace.vectors.c4,electrodePlace.startPoints.c4]=drillSpin(electrodePlace.angles.c4,cz,electrodePlace.c4,dicom);
[electrodePlace.vectors.t3,electrodePlace.startPoints.t3]=drillSpin(electrodePlace.angles.t3,cz,electrodePlace.t3,dicom);
[electrodePlace.vectors.t4,electrodePlace.startPoints.t4]=drillSpin(electrodePlace.angles.t4,cz,electrodePlace.t4,dicom);
[electrodePlace.vectors.oz,electrodePlace.startPoints.oz]=drillSpin(electrodePlace.angles.oz,cz,electrodePlace.oz,dicom);
[electrodePlace.vectors.pz,electrodePlace.startPoints.pz]=drillSpin(electrodePlace.angles.pz,cz,electrodePlace.pz,dicom);
[electrodePlace.vectors.fz,electrodePlace.startPoints.fz]=drillSpin(electrodePlace.angles.fz,cz,electrodePlace.fz,dicom);
[electrodePlace.vectors.fpz,electrodePlace.startPoints.fpz]=drillSpin(electrodePlace.angles.fpz,cz,electrodePlace.fpz,dicom);

[electrodePlace.vectors.f3,electrodePlace.startPoints.f3]=drillSpin(electrodePlace.angles.f3,cz,electrodePlace.f3,dicom);
[electrodePlace.vectors.f4,electrodePlace.startPoints.f4]=drillSpin(electrodePlace.angles.f4,cz,electrodePlace.f4,dicom);
[electrodePlace.vectors.p3,electrodePlace.startPoints.p3]=drillSpin(electrodePlace.angles.p3,cz,electrodePlace.p3,dicom);
[electrodePlace.vectors.p4,electrodePlace.startPoints.p4]=drillSpin(electrodePlace.angles.p4,cz,electrodePlace.p4,dicom);

%% load transducer specs

%freq=210000; % freq in hz
%period=0.5; % in seconds
voxelRes=0.001; % resolution in voxels, in meters
cmRes=10*voxelRes; % resolution in MM
Imax=3; % safety threshold in W/cm2 
% based on freq at 1 MHz
a.Blood=0.2;
a.Bone=6.9;
a.Brain=0.6;

%% calculate power spread
% start with direct lines

[electrodePlace.power.cz]=powerDis(electrodePlace.vectors.cz,a,Imax,cmRes);
[electrodePlace.power.c3]=powerDis(electrodePlace.vectors.c3,a,Imax,cmRes);
[electrodePlace.power.c4]=powerDis(electrodePlace.vectors.c4,a,Imax,cmRes);
[electrodePlace.power.t3]=powerDis(electrodePlace.vectors.t3,a,Imax,cmRes);
[electrodePlace.power.t4]=powerDis(electrodePlace.vectors.t4,a,Imax,cmRes);
[electrodePlace.power.oz]=powerDis(electrodePlace.vectors.oz,a,Imax,cmRes);
[electrodePlace.power.pz]=powerDis(electrodePlace.vectors.pz,a,Imax,cmRes);
[electrodePlace.power.fz]=powerDis(electrodePlace.vectors.fz,a,Imax,cmRes);
[electrodePlace.power.fpz]=powerDis(electrodePlace.vectors.fpz,a,Imax,cmRes);

[electrodePlace.power.f3]=powerDis(electrodePlace.vectors.f3,a,Imax,cmRes);
[electrodePlace.power.f4]=powerDis(electrodePlace.vectors.f4,a,Imax,cmRes);
[electrodePlace.power.p3]=powerDis(electrodePlace.vectors.p3,a,Imax,cmRes);
[electrodePlace.power.p4]=powerDis(electrodePlace.vectors.p4,a,Imax,cmRes);

powerModel=zeros(size(dicom));


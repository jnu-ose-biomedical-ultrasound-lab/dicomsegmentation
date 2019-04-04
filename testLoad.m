%--------------------------------------------------------------------------
 % testLoad.m

 % Last updated: April 2019, John LaRocco
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Loads a set of DICOM files and presents 3D results from it.
 % Currently hardcoded to use one dataset, but this can be changed.

%--------------------------------------------------------------------------


%% remove remnants of earlier iterations
close all;
clear;
clc;
% matlab dicom params: dinfo = dicominfo('TheFileName.dcm');
%spacing = dinfo.PixelSpacing;
%per_pixel_area = spacing(1) * spacing(2);
%num_white_pixels = nnz(binarized_img);
%total_white_area = num_white_pixels * per_pixel_area;


%% alternative: load fused dicom model from python

load('rawdicom.mat')
dicom=array;

%load('rawdicom.mat')
%dicom=array;


%% load dicom files using matlab
lbnd=417;
ubnd=622;
lbnd=0;
ubnd=284;

% i=1;
% file_name = 'I0000';
% file_name='image';
% range1=lbnd:1:ubnd;
% endstring=strcat(num2str(range1(i)),'.dcm');
% filename=strcat(file_name,endstring);
% dinfo = dicominfo(filename);
% spacing = dinfo.PixelSpacing;
% per_pixel_area = spacing(1) * spacing(2);
% X = dicomread(filename);
% firstFile=double(X);
% 
% dicomArray=zeros(size(firstFile,1),size(firstFile,2),length(range1));
% dicomArray(:,:,i)=squeeze(firstFile);
% 
% for i=2:(length(range1)-1);
% fileNo=range1(i);
% endstring=strcat(num2str(fileNo),'.dcm');
% filename=strcat(file_name,endstring);
% 
% X = dicomread(filename);
% ff=double(X);
% 
% dicomArray(:,:,i)=squeeze(ff);   
%     
%     
% end
% dicom=dicomArray;
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
%voxelRes=spacing; 
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


%% initialize power model
powerModel=zeros(size(dicom));

dicomModel=dicom;
powerModel(electrodePlace.cz(1),electrodePlace.cz(2),electrodePlace.startPoints.cz:(length(electrodePlace.power.cz)+electrodePlace.startPoints.cz-1))=(electrodePlace.power.cz);
powerModel(electrodePlace.c3(1),electrodePlace.c3(2),electrodePlace.startPoints.c3:(length(electrodePlace.power.c3)+electrodePlace.startPoints.c3-1))=(electrodePlace.power.c3);
powerModel(electrodePlace.c4(1),electrodePlace.c4(2),electrodePlace.startPoints.c4:(length(electrodePlace.power.c4)+electrodePlace.startPoints.c4-1))=(electrodePlace.power.c4);
powerModel(electrodePlace.t3(1),electrodePlace.t3(2),electrodePlace.startPoints.t3:(length(electrodePlace.power.t3)+electrodePlace.startPoints.t3-1))=(electrodePlace.power.t3);
powerModel(electrodePlace.t4(1),electrodePlace.t4(2),electrodePlace.startPoints.t4:(length(electrodePlace.power.t4)+electrodePlace.startPoints.t4-1))=(electrodePlace.power.t4);
powerModel(electrodePlace.oz(1),electrodePlace.oz(2),electrodePlace.startPoints.oz:(length(electrodePlace.power.oz)+electrodePlace.startPoints.oz-1))=(electrodePlace.power.oz);
powerModel(electrodePlace.pz(1),electrodePlace.pz(2),electrodePlace.startPoints.pz:(length(electrodePlace.power.pz)+electrodePlace.startPoints.pz-1))=(electrodePlace.power.pz);
powerModel(electrodePlace.fz(1),electrodePlace.fz(2),electrodePlace.startPoints.fz:(length(electrodePlace.power.fz)+electrodePlace.startPoints.fz-1))=(electrodePlace.power.fz);
powerModel(electrodePlace.fpz(1),electrodePlace.fpz(2),electrodePlace.startPoints.fpz:(length(electrodePlace.power.fpz)+electrodePlace.startPoints.fpz-1))=(electrodePlace.power.fpz);
powerModel(electrodePlace.f3(1),electrodePlace.f3(2),electrodePlace.startPoints.f3:(length(electrodePlace.power.f3)+electrodePlace.startPoints.f3-1))=(electrodePlace.power.f3);
powerModel(electrodePlace.f4(1),electrodePlace.f4(2),electrodePlace.startPoints.f4:(length(electrodePlace.power.f4)+electrodePlace.startPoints.f4-1))=(electrodePlace.power.f4);
powerModel(electrodePlace.p3(1),electrodePlace.p3(2),electrodePlace.startPoints.p3:(length(electrodePlace.power.p3)+electrodePlace.startPoints.p3-1))=(electrodePlace.power.p3);
powerModel(electrodePlace.p4(1),electrodePlace.p4(2),electrodePlace.startPoints.p4:(length(electrodePlace.power.p4)+electrodePlace.startPoints.p4-1))=(electrodePlace.power.p4);

if size(powerModel,3) > size(dicom,3),
powerModel(:,:,size(dicom,3))=squeeze(sum(powerModel(:,:,size(dicom,3):size(powerModel,3)),3));
powerModel=powerModel(:,:,1:size(dicom,3));
end


targetModel=powerModel;
targetModel(powerModel>0)=1;



%% save and export
save('electrodeInfo.mat','electrodePlace');
save('powerModel.mat','powerModel');
save('targetModel.mat','targetModel');

%% 2d plot
img=sum(targetModel,3);
figure();
imshow(img)
%% 3d plot

%% plot power
figure();
fv = isosurface(powerModel,.5);
p1 = patch(fv,'FaceColor','blue','EdgeColor','none');
view(3)

%% plot targets atop power
hold on;
fv = isosurface(targetModel,.5);
p1 = patch(fv,'FaceColor','red','EdgeColor','none');
view(3)
hold off;


% Slice function
% [X,Y,Z] = meshgrid(1:1:size(dicom,2));
% Z=Z(1:size(dicom,3));
% V=dicom;
% xslice = [];   
% yslice = [];
% zslice = 2:2:size(dicom,3);
% slice(V,xslice,yslice,zslice)

%% plot targets and power above shape
figure();
[x,y,z] = meshgrid(1:size(dicom,1),1:size(dicom,2),1:size(dicom,3));
data = dicom;
cdata = smooth3(rand(size(data)),'box',7);
p = patch(isosurface(x,y,z,data,10));
isonormals(x,y,z,data,p)
isocolors(x,y,z,cdata,p)
p.FaceColor = 'interp';
p.EdgeColor = 'none';
view(150,30)
daspect([1 1 1])
axis tight
camlight
lighting gouraud

hold on;
fv = isosurface(powerModel,.5);
p1 = patch(fv,'FaceColor','blue','EdgeColor','none');
view(3)

fv = isosurface(targetModel,.5);
p1 = patch(fv,'FaceColor','red','EdgeColor','none');
view(3)
hold off;

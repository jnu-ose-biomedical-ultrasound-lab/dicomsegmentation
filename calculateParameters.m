function electrodePlace=calculateParameters(dicom)
%--------------------------------------------------------------------------
 % calculateParameters.m

 % Last updated: April 2019, John LaRocco
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Calculates parameters of 3D DICOM data.
 % It is assumed that the dimensions are (X, Y, Z) or (Y, X, Z).

 % Inputs: 
 % dicom: A 3D matrix containing values from DICOM.
 
 % Outputs:
 % electrodePlace: A struct containing basic parameters, such as midpoints of non-zero X and y values. 
 
%--------------------------------------------------------------------------

%% Edge detection 
img=sum(dicom,3);
X=squeeze(sum(img));
Y=squeeze(sum(img'));

xnonzero=find(X~=0);
ynonzero=find(Y~=0);
xl=min(xnonzero);
xu=max(xnonzero);
xwidth=abs(xu-xl)+1;
%% calculate midpoints
yl=min(ynonzero);
yu=max(ynonzero);
ywidth=abs(yu-yl)+1;
midx=ceil(.5*xwidth);
midy=ceil(.5*ywidth);
cz=[midx,midy];
electrodePlace=[];
% note: This assumes the midpoint of all non-zero values is CZ. 
electrodePlace.cz=cz;

%% calculate ratios of other 10-20 international electrode distances
originVector=squeeze(dicom(midx,midy,:));
vXY=[xwidth,ywidth];
xLinePoints=[.1, .3, .5, .7, .9];
xZLines=[ceil(xwidth*xLinePoints); midy*ones(1,length(xLinePoints))]';
yLinePoints=[.1, .3, .5, .7, .9];
yZLines=[midx*ones(1,length(xLinePoints)); ceil(ywidth*yLinePoints)]';
quadrents=[.25, .3; .75, .3; .25, .7; .75, .7];
fused=[];
fused(1,:)=vXY(1).*quadrents(:,1);
fused(2,:)=vXY(2).*quadrents(:,2);
fused=ceil(fused');

%% add those values to our current object
electrodePlace.fused=fused;
electrodePlace.xZLines=xZLines;
electrodePlace.yZLines=yZLines;
electrodePlace.img=img;
electrodePlace.midx=midx;
electrodePlace.midy=midy;
electrodePlace.originVector=originVector;
end
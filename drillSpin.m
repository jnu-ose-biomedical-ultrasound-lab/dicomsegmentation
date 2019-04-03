function [vector,top]=drillSpin(theta,a,b,dicom)


%--------------------------------------------------------------------------
 % drillSpin.m

 % Last updated: April 2019, John LaRocco
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Return vector for given angle and data. 


 % Inputs: 
 % theta: Angle to rotate, an integer in degrees. Example: theta=33
 % a: A 2D array, both positive integers corresponding to one coordinate to compare between. First value is x, second is y. Example: a=[11,44]
 % b: The second 2D array, both positive integers corresponding to one coordinate to compare between. First value is x, second is why. Example: b=[4,6]
 % dicom: A 3D matrix of the slice to rotate. It is assumed the third dimension is Z. Example: dicom=array.
 
 % Outputs:
 % vector: 1D vector from the rotated matrix. 
 % top: A positive integer, corresponding to the index of the first non-zero value in the matrix. If none, defaults to top. 
 
%--------------------------------------------------------------------------


%% sort each 
x=sort([a(1),b(1)],'ascend');
y=sort([a(2),b(2)],'ascend');
dicomSlice=squeeze(dicom(x(1):x(2),y(1):y(2),:));
N = ndims(dicomSlice);
%% for 2D case
if N==2,
vector = radon(dicomSlice,theta);
topVec=find(vector~=0);

%% for 3d case
else
dicomSlice1=squeeze(dicom(x(1):x(2),a(2),:));
dicomSlice2=squeeze(dicom(a(1),y(1):y(2),:));
vector1 = radon(dicomSlice1,theta);
vector2 = radon(dicomSlice2,theta);

vector=zeros(max([length(vector1),length(vector2)]),1);

for i=1:length(vector1);
vector(i)=vector(i)+vector1(i);
end

for i=1:length(vector2);
if vector2(i) > vector(i),    
    
vector(i)=vector2(i);
end
end


end


%% find the first non-zero value
topVec=find(vector~=0);
if isempty(topVec);
top=1;    
else
top=topVec(1);    
end


end

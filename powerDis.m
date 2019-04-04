function [powerVec]=powerDis(vector,a,Imax,cmRes)

%--------------------------------------------------------------------------
 % powerDis.m

 % Last updated: April 2019, John LaRocco
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Model power dissipation through tissue. 

 % Inputs: 
 % vector: 1D vector, each corresponding to a line of voxels. 
 % a: Struct for attenuation values for bone, blood, and brain. 
 % Imax: Max intensity, in watts per cm^2
 % cmRes: a real positive value of how many cm each voxel corresponds to.
 
 % Outputs:
 % powerVec: 1D vector of same size as vector, showing power attentuation in W/cm2 at each voxel. 
 
%--------------------------------------------------------------------------

%% initialize
powerVec=zeros(1,length(vector));

for i=1:length(powerVec);
%Power=P times e^- (x*alpha)

%% tissue thresholds: soft tissue (0-12), unclotted blood (13-50), clotted blood (50-75), all blood (13-75), other tissue (75-299), bone (>300)
%vector(i)=round(vector(i));

R=round(vector(i));

if R <= 0 && R < 13,
a0=a.Brain;    
end

if R <= 13 && R < 76,
a0=a.Blood;    
end

if R <= 76 && R < 300,
a0=a.Brain;    
end

if R >= 300,
a0=a.Bone;    
end

eTerm=cmRes*a0*i;
coef=exp(-eTerm);

powerVec(i)=Imax*coef;
end

end

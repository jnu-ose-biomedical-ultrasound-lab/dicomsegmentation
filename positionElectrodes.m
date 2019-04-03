function electrodePlace=positionElectrodes(electrodePlace,F,V)
%--------------------------------------------------------------------------
 % positionElectrodes.m

 % Last updated: March 2019, John LaRocco
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: This file positions standard 10-20 international system
 % electrodes onto the DICOM-derived model. It can be adjusted for
 % different 'directions' 
 % If f=0, front is oriented towards lower index values.
 % If f=1, front is oriented towards higher index values.
 % If v=0, X axis is longer.
 % If v=1, Y axis is longer.
 
 % Inputs: 
 % electrodePlace: A struct with electrode positions and dimension
 % ratios.
 % F: A binary integer (0 or 1). Used for manual override.
 % V: A binary integery (0 or 1). Used for manual override. 
 
 % Outputs:
 % electrodePlace: Common electrode positions from 10-20 system added to the main struct. 
 

%--------------------------------------------------------------------------

%% determine the longer axis (x or y). Then determine front or back. 
img=electrodePlace.img;
midx=electrodePlace.midx;
midy=electrodePlace.midy;

% below means that Y is longer
if midy > midx;
    v=1;
    base1=img(1:midx,:); base2=img((1+midx):end,:);
   
% below means that X is longer     
else
     v=0;
     base1=img(:,1:midy); base2=img(:,(1+midy):end);
end



% Higher sum is assumed to be front. 
b1=sum(base1(:));  
b2=sum(base2(:));

if b1 >= b2; % front is oriented towards higher index values
f=1;
else % front is oriented towards lower index values
f=0;
end
%% Manual override
if nargin > 1;

if isempty(V);

else
v=V;    
end


if isempty(F);

else
f=F;    
end

end




%% logic tree for coordinate fitting

xZLines=electrodePlace.xZLines;
yZLines=electrodePlace.yZLines;
fused=electrodePlace.fused;
electrodePlace.bearing=[f,v];

%% facing north
if f==1 && v==1; % front is near higher indexes, and Y is longer axis
electrodePlace.c3=xZLines(2,:);
electrodePlace.c4=xZLines(4,:);
electrodePlace.t3=xZLines(1,:);
electrodePlace.t4=xZLines(5,:);
electrodePlace.oz=yZLines(1,:);
electrodePlace.pz=yZLines(2,:);
electrodePlace.fz=yZLines(4,:);
electrodePlace.fpz=yZLines(5,:);

electrodePlace.f3=fused(3,:);
electrodePlace.f4=fused(4,:);

electrodePlace.p3=fused(1,:);
electrodePlace.p4=fused(2,:);
electrodePlace.facing='north';
end
%% facing south
if f==0 && v==1; % front is near lower indexes, and Y is longer axis
electrodePlace.c3=xZLines(4,:);
electrodePlace.c4=xZLines(2,:);
electrodePlace.t3=xZLines(5,:);
electrodePlace.t4=xZLines(1,:);
electrodePlace.oz=yZLines(5,:);
electrodePlace.pz=yZLines(4,:);
electrodePlace.fz=yZLines(2,:);
electrodePlace.fpz=yZLines(1,:);
electrodePlace.f3=fused(2,:);
electrodePlace.f4=fused(1,:);
electrodePlace.p3=fused(4,:);
electrodePlace.p4=fused(3,:);
electrodePlace.facing='south';
end

%% facing west
if f==1 && v==0; % front is near higher indexes, and X is longer axis
electrodePlace.c3=yZLines(2,:);
electrodePlace.c4=yZLines(4,:);
electrodePlace.t3=yZLines(1,:);
electrodePlace.t4=yZLines(5,:);
electrodePlace.oz=xZLines(5,:);
electrodePlace.pz=xZLines(4,:);
electrodePlace.fz=xZLines(2,:);
electrodePlace.fpz=xZLines(1,:);
electrodePlace.f3=fused(1,:);
electrodePlace.f4=fused(3,:);
electrodePlace.p3=fused(2,:);
electrodePlace.p4=fused(4,:);
electrodePlace.facing='west';
end
%% facing east
if f==0 && v==0; % front is near lower indexes, and X is longer axis
electrodePlace.c3=yZLines(4,:);
electrodePlace.c4=yZLines(2,:);
electrodePlace.t3=yZLines(5,:);
electrodePlace.t4=yZLines(1,:);
electrodePlace.oz=xZLines(1,:);
electrodePlace.pz=xZLines(2,:);
electrodePlace.fz=xZLines(4,:);
electrodePlace.fpz=xZLines(5,:);
electrodePlace.f3=fused(4,:);
electrodePlace.f4=fused(2,:);
electrodePlace.p3=fused(3,:);
electrodePlace.p4=fused(1,:);
electrodePlace.facing='east';
end

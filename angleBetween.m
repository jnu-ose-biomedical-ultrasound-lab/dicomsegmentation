function theta=angleBetween(cord1,cord2)


%--------------------------------------------------------------------------
 % angleBetween.m

 % Last updated: April 2019, John LaRocco
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Calculates angle coordinates for different electrode pairs, based on electrode placement.


 % Inputs: 
 % cord1: An array with two integer coordinates, representing a top-down view on
 % the skull. Example: cord1=[45,32] 
 
 % cord2: An array with two integer coordinates, representing a top-down view on
 % the skull. Example: cord2=[65,22] 
 
 % Outputs:
 % theta: The angle between cord1 and cord2, expressed in degrees. 
 
%--------------------------------------------------------------------------


%% separate each 

x1=cord1(1);
y1=cord1(2);

x2=cord2(1);
y2=cord2(2);

theta = rad2deg(atan2(x1*y2-x2*y1,x1*x2+y1*y2));



end

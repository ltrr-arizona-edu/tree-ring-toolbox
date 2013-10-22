function R=dms2rad(C)
% Angle in deg, min, seconds to angle in radians

%*********************  INPUT
%
% C (mC x 3)  angle in degrees, minutes, seconds


%************  OUTPUT
%
% R (mR x 1)  equivalent angle in radians



R=(C(:,1) + C(:,2)/60 + C(:,3)/3600) * pi/180;
% latchg.m  latitude-longitude change convention

% Given a file of longs and lats in degrees and minutes, gives
% an array Z of same in decimal degrees.  Assumes input file has
% minus sign before west longitude degrees.  (Not before the minutes).
% x-coord in Z will be negative to the west for convention of 
% interpreting plots.

load c:\projs\aa3\nsf.lat
z=nsf;
z(:,2)=z(:,2) ./ 60;
x=z(:,1) +z(:,2);


z(:,4)=z(:,4) ./ 60.0;

y=z(:,3)- z(:,4);

z=[x y];


function [CEQ,CE2,CE3,CS1,CS2,CS3]=stat1a(IX,EX,yrs,nms)
%
% DESCRIPTION : stats1
% Mean, variance and mean sensitivity of core indices 
%
% INPUTS  :  NONE
%
% OUTPUTS :  NONE
%
% USER WRITTEN FUNCTIONS NEEDED :
% 	- BASIC1.M
% 	- RTREE1.M
%	- RTREE2.M
%	- MEANSEN1.M
%___________________________________________________________



% Get the mean and the std dev for IX
[mn,std,n] = basic1(yrs,IX);
ms = meansen1(IX,yrs);

% Build the CS1 matrix
CS1=[(1:length(mn))' yrs(:,1:2) n mn ms std];

% Get the mean and the std dev for EX
[mn,std,n] = basic1(yrs,EX);
ms = meansen1(EX,yrs);

% Build the CE1 matrix
CE1=[(1:length(mn))' yrs(:,1:2) n mn ms std];

% Get the Correlation coefficients
CS2 = rtree1(IX,nms,yrs);
CE2 = rtree1(EX,nms,yrs);

% Get the average correlation coefficients
CS3 = rtree2(CS2);

% Get the average correlation coefficients
CE3 = rtree2(CE2);



% End of file

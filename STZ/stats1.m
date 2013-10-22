function [CE1,CE2,CE3,CS1,CS2,CS3]=stats1(IX,EX,yrs,nms)
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

%*** UW FUNCTIONS CALLED
% 
% meansen1
% basic1
% rtree1
% rtree2
% maskind
% treenum

[ns,dum]=size(nms);

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

% Get the Correlation coefficients between and within trees
CS2 = rtree1(IX,nms,yrs);
CE2 = rtree1(EX,nms,yrs);

% Get the average correlation coefficients
CS3 = rtree2(CS2); % standard index
CE3 = rtree2(CE2); % residual index


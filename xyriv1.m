function [lat,long]=xyriv1(DNline)
% xyriv1:  get line file for Sacramento River, California, from worldlo
% CALL: xyriv1;
%
% Meko 4-3-98
%
%****************
%
L1 = DNline(1).long >=-130 & DNline(1).long<=-119;
L2 = DNline(1).lat >=35 & DNline(1).lat<=45;
L3=L1 & L2;

lata = DNline(1).lat(L3);
longa = DNline(1).long(L3);

disp('here');


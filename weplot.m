% weplot.m    west to east plotting of tree-climate median correlations

% Scatterplot of r tree vs climate against either long 
% May need to modify, depending on (1) season, (2) max or medn corr, 
%   (3) period for anlysis--1 fewer year for autumn.

% Possible change lines flagged with "hhh".  These lines would change
% most often.  Other lines could change if you use search distances other
% than 110 km

%  copy this file over to your working directory before running.  That
% way avoid adding a lot of silly paths.



% preload: nsf.lat,  pwi110.dat,  distp110

treedeg=nsf;

R1=pfa110;    % hhh


axis([-125  -60  -.7 .7]);  % hhh  [25 50] or [-125  -60]

x=treedeg(I2,2);    %  Only those sites with corrln coefs in MR...
y=treedeg(I2,1);


z=R1(:,4);   %  4 or 5, depending on if median or max


plot(x,z,'*',[-125 -60],[0 0],'--r') % hhh [25 50] or [-125 -60]

text(-124,-.55,'SEP-NOV');    % hhh  ANY OF 4 SEASONS
text(-124,-.65,'SEARCH RADIUS 110 KM'); % [26 -.55] or [124 -.55] 

xlabel('LONGITUDE (-DEG)') 
ylabel('MEDIAN R');   % hhh  MEDIAN OR MAXIMUM
title('TREE RING-CLIMATE CORRELATION, 1906-1979')  ;   % hhh  1906 OR 07

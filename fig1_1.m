% fig1_1.m  

% for time series plots of factor scores based on original data
% and on series variance-stabilized and whitened

%  To run for first time in session, uncomment the following block.
%  For subsequent runs, comment it out again.  Change columns in the plot
%  command, and the title info to match whatever factors you want plotted.
%  Since the year in res and org is col 1, factor 1 is col2,
%    factor 2 is col 3, etc.

%  Factor 4 in res most closely resembles factor 3 in org, and vice
%  versa.   Other factors correspond 1 to 1.


%########### comment-out block ######################

% load c:\projs\aa3\orgscore.dat
% load c:\projs\aa3\resscore.dat

%res=resscore;
%org=orgscore;
%yr=orgscore(:,1);

%clear rescore orgscore

%############# end of comment-out block #############

clg

plot(yr,org(:,3), yr,res(:,3));
title('FACTOR 2')

pause

plot(org(:,3),res(:,3))

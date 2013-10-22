function Z =  plotmk(names3,x3,yrs);

% Double mass plots and Mann Kendall statistics for 3 station-pairs
%
% D. Meko     11-23-93
%
%
% names3 (3 x 12)  string matrix of station names and season
% x3     (mX3 x 3) seasonal or monthly ppt at 3 stations
%
% Z  (3 x 2)  col 1 is tau value, col 2 is 95% prob point
%		row 1 is station 1 vs station 2
%		row 2 is station 1 vs station 3
%		row 3 is station 2 vs station 3
%
% yrs (2 x 1) starting and ending years 
%
%
%*********  FUNCTIONS NEEDED --- mannkend.m
%
%
%**********   USE   *****************************************
%
% Assume have pulled relevant station time series, same length, into
% the 3-col matrix x3.  Assume have pulled station/season names into
% 12-col string matrix names3.   For example, names3 might be:
%
%  tucson nv-ap
%  phoenx nv-ap
%  bisbee nv-ap
%
% and x3 would hold the nov-april ppt totals, 1920-1990 at each
% station.
%
% Run the function.  The returned values in Z give the Mann-Kendall
% tau and 95% prob points. If tau is greater than the 95% point, 
% one series in inhomogeneous relative to the other at the 0.05
% confidence level.  The function also plots 4 frames in the graphics
% window.  Three frames are cumulative ppt vs cumulative ppt for each
% station pair.  The lower right frame has unlabeled time series plots
% for the three precip series.  The is a quick and dirty plot to scan
% for obvious evidence of change in level of time series.  Yellow is
% series 1,  cayan is series 2, purple is series 3.  Colors may vary
% depending on system.

Z=zeros(3,2);  % preallocate

I1=[1 2; 1 3; 2 3];  % specify sequence of testing of pairs
mnp=[2 2 1;  2 2 2;  2 2 3];  % specify plot frames for subplot

% Compute cumulative sums of ppt at each station
ys=[cumsum(x3(:,1))  cumsum(x3(:,2)) cumsum(x3(:,3))]; 

% Loop thru each pair of stations
for k = 1:3;
	i=I1(k,1);
	j=I1(k,2);
	x=x3(:,i);  % x and y are the original series, in inches or whatever
	y=x3(:,j);
	[tau,p95]=mannkend(x,y);  % compute mann kendall stat
	Z(k,:)=[tau p95];

	subplot(mnp(k,1),mnp(k,2),mnp(k,3));  % specify one of the four frames
	plot(ys(:,i),ys(:,j),'*');  % plot cum sum of one series vs anoth
	xlabel(names3(i,:));  % label for x-axis of cum-sum plot
	ylabel(names3(j,:));  % label for y-axis

	% Compute some minima and maxima needed to locate text on plot
	mnx2=min(ys(:,i));
	mxx2=max(ys(:,i));
	mny2=min(ys(:,j));
	mxy2=max(ys(:,j));

	% Compute plotting points for text
	xp1=mnx2+(mxx2-mnx2)/10;
	yp1=mny2+9.5*(mxy2-mny2)/10;
	text(xp1,yp1,['tau = ',num2str(tau)]);
	xp2=mnx2+3*(mxx2-mnx2)/10;
	yp2=mny2+1*(mxy2-mny2)/10;
	text(xp2,yp2,['p95 = ',num2str(p95)]);
end
	subplot(2,2,4);  % specify frame for unlabled time series plots
	yr=(yrs(1):yrs(2))';   % year vector
	plot(yr,x3(:,1),yr,x3(:,2),yr,x3(:,3));

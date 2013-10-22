function [p25,nl,y]=fltplay1(x,yr)
%[p25,nl,y]=fltplay1(x,yr)
% Trial and error fits of FIR filters.  
% D. Meko 11-5-92
%*******  INPUT ARGS
%
% x time series, col vect
% yr corresp years
%**********  OUTPUT ARGS
%
% p25  (1x1) period (yr) of 25% freq response
% nl  (1x1)  length of filter
%************  COMMENTS **********************
%
% Usually run before firdmm1 to arrive at settings for p25 and nl
%


p25=input('25% freq response (yr):  ');
nl= input('Desired number of weights in filter: ');
nex=nl; % time series will be appended with this many years of 
%      median values on each end before passing thru filter


n=nl-1;  
Wn=2.0 / p25;  % Cutoff frequency on scale such that 0 is 
	% 0 per year, and Wn=1 is 0.5 per year (Nyquist).  At
	% this frequency, response of filter b1 is 0.5, and
	% response of b=conv(b1,b1) is 0.25.

b=fir1(n,Wn);  % Note: before doubling

xm=median(x);
x1=[xm(ones(nex,1),1); x; xm(ones(nex,1),1)]; % extend both ends nex yr
%	with medians

y1=filtfilt(b,1,x1);
y=y1(nex+1:length(y1)-nex);
plot(yr,x,yr,y);
title(['25% response = ',num2str(p25),'  length= ',int2str(nl)])
pause

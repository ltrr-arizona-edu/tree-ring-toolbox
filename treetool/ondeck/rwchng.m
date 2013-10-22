function y=rwchng(x)
% ring-width change series for use with xdate.m
% Steps 
%  1. 9-wt hamming low-pass filter ring width
%  2. departures from low-pass series curve (mm)
%  3. first difference of departures (mm) to accentuate change in rw
%  4. ratio of first diff of departures to low-pass series, to
%	scale to local mean of rw


%******   INPUT ARGS
%
% x (mx x 1)  ring-width series


%*****   OUTPUT ARGS
%
% y (my x 1)  transformed ring width
%     my=mx-1


mx=length(x);
% Design low-pass Hamming Filter
b=fir1(8,.1);

% Filter the ring width
g=filtfilt(b,1,x);  


% Departures from smooth curves
d=x-g;

% First differences of departures
f=diff(d);

% Standardized to value of smooth curve in previous year
gs=g(1:mx-1);
y=f ./ gs;
t=(1:length(y))';
%plot(t,x(2:length(x)),t,y,t,gs)';

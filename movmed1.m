function [y,yry]=movmed1(x,yr,n1,kopt)
%
% Compute moving median of time series 
%
% D Meko 5-1-97
%
% Built on movsum1.m template. 
%  
%*******  INPUT ARGS
%
% x time series, col vect
% yr corresp years
% n1   number of observations in window
% kopt (1 x 1)i options
%     kopt(1):  ==1 plot year at end of period
%				==2 plot year at middle of period
% 
%**********  OUTPUT ARGS
%
% y   (?x1) moving median
% yry (? x1) plotting year vector for y 
%		Either at end or middle of n1-year periods, depending on
%		kopt(1)


% Compute length of filtered series
ny=length(yr)-n1+1;

% Compute 'plotting year' vector
yrinc=(0:ny-1)'; % will add this to start year for year vector
if kopt(1)==1; % plot value at ending year
	yr1=yr(1)  + n1 -1;
	txt1='Ending Year';
elseif kopt(1)==2; % plot at central point in n1-year period
	yr1=yr(1)+(n1-1)/2;
	txt1='Middle Year'
else
	error('kopt(1) must ==1 or 2');
end
yry=yr1+yrinc; % vector of plotting years

if n1==1;
   y=x;
   return
end


mx=length(x);
t=(1:mx)'; % time vector for original series

% Compute how many windowed periods
npers = mx-n1+1;

% Compute rv of start index of each period
i1=[1:(mx-n1+1)];

% Compute matrix of indices to each windowed period -- each col represents a period
j1=(0:(n1-1))'; % cv
J1=repmat(j1,1,npers);
I1=repmat(i1,n1,1); 
I2=I1+J1;

% Compute the medians
y=(median(x(I2)))';

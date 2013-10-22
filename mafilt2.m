function [y,yry]=mafilt2(x,yr)

% Compute evenly weighted moving average of time series.  

% D. Meko 2-26-93


%*******  INPUT ARGS
%
% x time series, col vect
% yr corresp years
% m   number of weights

%**********  OUTPUT ARGS
%
% y   (?x1) smoothed time series
% yry (? x1) years for y


nl= input('Desired length of moving average: ');

% Compute weights

weight=1/nl;
b=weight(:,ones(nl,1));

y1=filter(b,1,x);  % filtered series, before shifting and truncating

% Compute length 

ny=length(yr)-length(b)+1;
yrinc=(0:ny-1)';
yr1=yr(1)+(length(b)-1)/2;
yry=yr1+yrinc;

y=y1;
y(1:length(b)-1)=[];

plot(yr,x,yry,y);
title(['Moving average of length ',int2str(nl)])
pause

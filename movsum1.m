function [y,yry]=movsum1(x,yr,n1,kopt)
%
% Compute moving sum of time series 
%
% D Meko 7-26-96
%
% Built on mafilt2.m template.  Needed for time series of 
% frequency of wet and dry years in ten-year periods
%  
%*******  INPUT ARGS
%
% x time series, col vect
% yr corresp years
% n1   number of weights
% kopt (1 x 2)i options
%     kopt(1):  ==1 plot year at end of period
%				==2 plot year at middle of period
%		kopt(2)	==1 do not plot within this call
%				==2 plot within this call
% 
%**********  OUTPUT ARGS
%
% y   (?x1) moving sum
% yry (? x1) plotting year vector for y 
%		Either at end or middle of n1-year periods, depending on
%		kopt(1)


% Build unit kernel
b=ones(1,n1);

% filter series
y1=filter(b,1,x);  


% Compute length of filtered series
ny=length(yr)-length(b)+1;

% Compute 'plotting year' vector
yrinc=(0:ny-1)'; % will add this to start year for year vector
if kopt(1)==1; % plot value at ending year
	yr1=yr(1)  + length(b) -1;
	txt1='Ending Year';
elseif kopt(1)==2; % plot at central point in n1-year period
	yr1=yr(1)+(length(b)-1)/2;
	txt1='Middle Year'
else
	error('kopt(1) must ==1 or 2');
end
yry=yr1+yrinc; % vector of plotting years

% Truncate start-up years of time series
y=y1;
y(1:length(b)-1)=[];

% Optionally plot original and smoothed series
if kopt(2)==1; % no plot
elseif kopt(2)==2; 
	figure
	subplot(2,1,1);
	plot(yr,x);
	title('Original Series')
	subplot(2,1,2);
	plot(yry,y);
	title(['Moving sum of length ',int2str(n1)])
	xlabel(txt1);
	kclose=input('Close the Window? Y/[N]:   ','s');
	kclose=upper(kclose);
	if isempty(kclose);	
		kclose='N';
	elseif kclose~='Y' & kclose~='N' 
		error('kclose must by Y or N')
	end
	if kclose=='Y'
		close all;
	else
		disp('Plots are in figure window')
	end
else
	error('kopt(2) must be 1 or 2')
end

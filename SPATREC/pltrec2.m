function [e]=pltrec2(yrc,yc,yr,y,pebar,modnum,hf)
% pltrec2(yr,y,pebar,modnum)
%
% tsp of long-term reconstruction, including recon for calib period.
% horiz line at the calib-period mean of recon data
% lines at +- rmse from horiz line.  rmse is from validation.
% asterisks at extrapolations
mn=mean(yc);
a=NaN;

yrtrunc='N'; % truncate beginning year of reconstruction plot?
scalef=1.0;  % scaling factor to make figure in cm for Holocene

tocm = input('Convert inches to cm for plot? Y/[N]  ','s');
if isempty(tocm) | tocm=='N' | tocm=='n',
	tocm='N';
elseif tocm=='Y' | tocm=='y',
	tocm='Y';
else 
	error('Y,y,N or n are only allowable responses')
end


yrtrunc=input('Truncate start year of plot? Y/[N] ','s');
if isempty(yrtrunc) | yrtrunc=='N' | yrtrunc=='n'
	yrtrunc='N';
elseif yrtrunc=='Y' | yrtrunc=='y',
	yrtrunc='Y';
else
	error('Y,y,N or n are only allowable responses')
end


	
if tocm=='Y',
	scalef=2.54;
	yc=yc*scalef;
	pebar=pebar*scalef;
	y=y*scalef;
	mn=mn*scalef;
end


% Compute time series of rmse
for n = 1:length(modnum);
	e(n) = pebar(modnum(n)); % root-mean-square prediction error
	e=e';
end

% Flags for extrapolations
f=max(e)+.05;
f=f(ones(length(yr),1),1);
f = f(hf);

% Prompt for start year of plot, then
% Calculate how many years to lop off start of each time series
if yrtrunc=='Y',
	yrgo=input('Start year for plot: ');
	if yrgo < yr(1) | yrgo >yr(length(yr)),
		error('Specified start year out of range')
	end
	nout = yrgo-yr(1); % number to drop
	% Drop the leading nout years from plotted series
	if nout > 0;
		yr(1:nout)=[];
		y(1:nout)=[];
		e(1:nout)=[];
		hf(1:nout)=[];
	end
end

yr3=yr(hf);

h1=plot(yr,y,yrc,yc,[yr(1) yrc(length(yrc))],[mn mn])
h2=line(yr,mn+e);
h3=line(yr,mn-e);
h4=line(yr3,y(hf))

set(h4,'LineStyle','*','MarkerSize',[4])

v1=axis;
xlabel('Year')
if tocm=='Y',
	ylabel('PPT (cm)')
else
	ylabel('PPT (in)')
end
title('RECONSTRUCTED COOL-SEASON PRECIPITATION')





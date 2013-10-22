function [e]=pltrec1(yrc,yc,yr,y,pebar,modnum,hf,ktrans)
% pltrec1(yr,y,pebar,modnum)
%
% tsp of long-term reconstruction, including recon for calib period.
% horiz line at the calib-period mean of recon data
% lines at +- rmse from horiz line.  rmse is from validation.
%
% See spatrec1.m for instructons

mn=mean(yc);
a=NaN;

% Compute time series of rmse
for n = 1:length(modnum);
	e(n) = pebar(modnum(n)); % root-mean-square prediction error
	e=e';
end


% Check current figure background color and change if necessary
curcolor=get(gcf,'Color');
black=all(curcolor==0);
if ktrans==0 & black==1; % OK for black bkgd slide
elseif ktrans==0 & black==0; % change background
	whitebg;
elseif ktrans==1 & black==0; % OK for white background transparency
elseif ktrans==1 & black==1; % change background
	whitebg
else
	error('unknown combination of ktrans and black');
end


% Plot figure
	h1=plot(yr,y,yrc,yc,[yr(1) yrc(length(yrc))],[mn mn])
	h2=line(yr,mn+e);
	h3=line(yr,mn-e);

% Change colors if needed
if ktrans==0; % want black background. no changes needed
elseif ktrans==1; % for transparency
	set(h1(1),'Linewidth',2,'Color','m')
	set(h1(2),'Linewidth',1,'Color','k')
	set(h1(3),'Color','k')
	set(h3,'Color','r')
else
	error('Invalid settind for ktrans')
end


v1=axis;
ylabel('PPT (in)')
title('RECONSTRUCTED COOL-SEASON PRECIPITATION')
legend('Reconstruction','Observed','Calib-Pd Mean','1 RMSE')




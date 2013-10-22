function yup = setylim2 (digs,d,youter)
% set Ylim to reserve upper part of plot for text (see pltext2.m)
% D Meko,  10-15-94
%
%
%************************  IN ARGS ******
%
% digs (1 x 1)i   number of digits to right of decimal point in 
%	the y-axis tic labels.  For example, if y axis will have
%	labels at [0.0  0.1  0.2  ...], digs=1.  Or, if labels
%	are at [0.12 0.13  0.14  ... ], digs = 2.
%	digs is needed to properly round the ylimits
% d (1 x 1)r   decimal fraction of original ylimits by which to expand
%	ylimits to reserve blank area for label info.  
%	Typicall, something like 0.1 (10%)
% youter (1 x 2)r  minimum and maximum of data values 
%
%
%************************************ OUT ARGS ********
%
% yup (1 x 1)r   y-value (value on y-axis) marking lower boundary
%	of text region.  yup is needed in subsequent call to pltext2.m
%
%
%************************ USER-WRITTEN FUNCTIONS NEEDED -- NONE
%
%
%**************************  NOTES ********************************
%
% Example call:
%
%   yup = setylim2 (1,.1,[min(y) max(y)])
% 
%   Says expand the ylimits by 10% (0.1) and consider the y-tick labels
%   to have 1 place to right of decimal point.    



ylo=youter(1);
yhi=youter(2);

% Calc Ylim values that would give round values for lowest and
% highest tick labels
ymin2 = (floor(ylo * (digs*10)))/(digs*10);
ymax2 = (ceil (yhi * (digs*10)))/(digs*10);

% Calc range of this y axis
dely = ymax2-ymin2;

% Calc upper y axis limit for expanded axis
ymin3 = ymin2;
ymax3 = ymax2 + dely*d;
yup=ymax2;

% Set the y-axis limits
set(gca,'Ylim',[ymin3,ymax3])
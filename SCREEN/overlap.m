function [x,y,yr]=overlap(v,w,YRS)
% [x,y,yr]=overlap(v,w,YRS)
%
% Given two time series, column vectors covering possibly different
% periods, get the correponding parts of the series for the common 
% period. Input vectors v and w can have leading or trailing 
% NaN's
%
% D Meko 4-27-96
%
%******************* IN ARGS ************************
%
% v (mv x 1)r one input time series 
% w (mw x 1)r the other series
% YRS (2 x 2)i  first and last years of the vectors
%		v (row 1) and
%		w (row 2)
% 
%
%********************* OUT ARGS ***********************
%
% x (mx x 1)r subset of v for the common period
% y (my x 1)r subset of w for the common period; mx==my
% yr (mx x1)i year vector for x and y
%
%***************** USER WRITTEN FUNCTINS NEEDED -- NONE
%
%***************** NOTES ***************************
%
% The non-NaN values in v and w must be on either end, or
% both ends, not internal -- i.e., not imbedded in the 
% sequence of valid data

yrv = (YRS(1,1):YRS(1,2))';
yrw = (YRS(2,1):YRS(2,2))';
L1 = ~isnan(v);
L2 = ~isnan(w);
if ~any (L1) | ~any(L2),
	error('v or w are all NaN');
end
v = v(L1);
w= w(L2);
yrv = yrv(L1);
yrw = yrw(L2);

if ~all(diff(yrv)==1) | ~all(diff(yrw)==1)
	error('v or w have internal NaN values')
end

% Compute start and end year of common period
yrgo = max([min(yrv) min(yrw)]);
yrsp = min([max(yrv) max(yrw)]);

% Check that some overlap between the two non-NaN parts of the 
% two series
if (yrsp-yrgo)<0,
	error ('No overlap in valid parts of the two series')
end

% Pointers to desired part of each series
Lv = yrv>= yrgo & yrv <=yrsp;
Lw = yrw>=yrgo & yrw<= yrsp;

yr = (yrgo:yrsp)'; % year vector for output
x = v(Lv);
y = w(Lw);



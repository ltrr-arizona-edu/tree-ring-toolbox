function [a,b,xends,tends,z1,z2,z3]=triline(x,yr1)
% [ao,b]=triline(x,yr1)
%
% A resistant line from three groups
% 
% D. Meko   3-8-96
%
%
%******************** IN ARGS ****************************
%
% x (mx x 1)r    time series
% yr1 (1 x 1)i   first year of x
%
%
%********************* OUT ARGS **********************
%
% a (1 x 1)r  intercept of line
% b (1 x 1)r  slope of line
%
%***************** USER WRITTEN FUNCTIONS CALLED ************
%
% trigroup.m  -- divide time series in three groups
%
%
%***************  NOTES **************************************
% SOURCE: Hoaglin, D.C, Mosteller, F., and Tukey, J.W., 1982, Understnding
% Robust and Exploratory Data Analysis, John Wiley & Sons, 130.
%
% Time series split into early, middle, late thirds.
% Median of x and time (year) computed for each third
% Line fit to those summary points
%
%**********************************************************

% checks
[mx,nx]=size(x);
if mx<2 | nx~=1,
	error('x must be a column vector')
end
[yr1m,yr1n]=size(yr1);
if yr1m~=1 | yr1n ~=1,
	error('yr1 must be a scalar')
end

yr = (yr1:yr1+mx-1)';  % year vector


% Get start years of segments
t1=trigroup(x,yr1);

% Compute ending years of segments
t2=[t1(2)-1   t1(3)-1   yr(length(yr))];

L1 = yr>=t1(1) & yr<=t2(1);
L2 = yr>=t1(2) & yr<=t2(2);
L3 = yr>=t1(3) & yr<=t2(3);

% Get segments of the time series and years
x1=x(L1);
k1=yr(L1);
x2=x(L2);
k2=yr(L2);
x3=x(L3);
k3=yr(L3);

% get medians
z1=median([k1 x1]);
z2=median([k2,x2]);
z3=median([k3,x3]);

xL=z1(1);  % These in Hoaglin notation, p. 132
yL=z1(2);
xM=z2(1);
yM=z2(2);
xR=z3(1);
yR=z3(2);

% Compute slope
b = (yR-yL)/(xR-xL);

%  Intercept
a = (1/3)* ((yL-b*xL)+(yM-b*xM)+(yR-b*xR));


% Compute x,y values for the fitted line
zgo = a + b * yr1;
zsp = a + b * yr(mx);
xends=[zgo zsp];
tends=[yr1  yr(mx)];

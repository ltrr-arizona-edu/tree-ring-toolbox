function [x,y]=getxy(X,tabx,taby)
% getxy.m  using lookup table, convert lat, long to x,y

% D Meko   4-94
%  Last revised 5-9-94: converted to function
% 
%
%*********** IN ARGS ***********
%
% X (mX x 3) longitude and latitude for mX points.  Col 1 is seq # for point;
%   col 2 is long in decimal deg;  col 3 is lat in decimal deg
% tabx, taby:  two lookup tables typically stored in a .mat file that goes along with the
%   .dwg and other dedicated drawing files.  For example, cf17.dwg, the 1:1,000,000 map of
%    the northern Plains has a file cf17.mat, which holds tabx and taby.  The structure of
%   tabx and taby  is as follows.  Col 1 is longitude (x-direction), negative and decreasing
%   in absolute value with row number.  Col 2 is latitude (y-direction), positive, and increas
%   ing with col number.  For example, here is a part of tabx:
% 
%   -105.0000   40.0000   41.0000   42.0000   43.0000   44.0000
%   -105.0000   13.7066   13.6066   13.8059   14.0403   14.2583
%   -104.0000   17.0074   16.9074   17.0587   17.2241   17.3837
%   -103.0000   20.3144   20.2144   20.3134   20.4165   20.5249
%   -102.0000   23.6258   23.5258   23.5743   23.6147   23.6807
%   -101.0000   26.9492   26.8492   26.8460   26.8174   26.8151
%   -100.0000   30.2220   30.1220   30.0585   29.9907   29.9306
%    -99.0000   33.5286   33.4286   33.3153   33.1955   33.0898 
%
%   Table entries for tabx are the x plotting points in the dwg drawing, in drawing units.
%   Table entries for taby are the y plotting points
%
%   Note that tabx(1,1) and tabx(2,1) are the same.  Col 1 must be
%   non-decreasing downward.  Row 1 must be non-decreasing to right. The convention
%   I use agrees with this rule.  Be sure that all the data points to be interpolated
%   fall fall within the lat range between col2 and the next-to
% 
%
%************ OUT ARGS *************
%
% x (mX x 1) -- x-coord of point in drawing units
% y (mX x 1) -- y-coord of point in drawing units
%
%
%
%**************  USER FUNCTIONS NEEDED -- none
%
%
%
%


[mx,nx]=size(X);
for i=1:mx;
	x(i)=table2(tabx,X(i,2),X(i,3));
	y(i)=table2(taby,X(i,2),X(i,3));
end
x=x';
y=y';

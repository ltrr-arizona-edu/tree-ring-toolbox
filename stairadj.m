function [X1,Y1]=stairadj(X,Y,baselev)
% [X1,Y1]=stairadj(X,Y)
% 
% Adjust "stairs" left-side output for correct time series plotting
%
% D Meko 11-27-96
%
%
%**************** IN ARGS 
%
% X,Y -- obtained by calling function using call to stairs. For example,
%	[X,Y]=stairs(yr,z)
%	X would be the abscissa values for stair plot, Y the ordinate
% baselev (1 x 1)r  base level of the stairs -- usually zero
%
%**************** OUT ARGS
%
% X1,Y1 -- adjusted x,y plotting points for X and Y.  
%
%
%************** NOTES
%
% Needed to "pretty up" start and end of stairs plots of time series.
%
% Typical use: Have timeseries x and years yr.  Run 
% [X,Y]=stairs(yr,x).  Call [X1,Y1=stairadj(X,Y,baselev).
% Then -- plot(X1,Y1),

% adjust stairs variables for correct positioning
xbeg=X(1);
ybeg=baselev;
xend=[X(length(X))+1;  X(length(X))+1];
yend=[Y(length(Y)); baselev];
X1=[xbeg; X; xend]-0.5;
Y1=[ybeg; Y; yend];


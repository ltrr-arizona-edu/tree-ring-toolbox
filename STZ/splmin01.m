function y=splmin01(f,p);
% splmin01: minimizing function to find period of 50% spline response given spline parameter
% y=splmin01(f,p);
% Last revised 9-9-99
%
% Called by spline50.m used to compute the 50% frequency response period of a 
% spline.  Useful typically if have specified your spline by some method other than
% specifying the period of its 50% frequency response
%
%*** INPUT
%
% f (1 x 1)r  frequency (1/yr) of starting point for iteration
% p (1 x 1)r  spline parameter
%
%*** OUTPUT 
%
% y (1 x 1)r  dummy minimizing value


A=6*(cos(2*pi*f)-1).^2;
B=(cos(2*pi*f) + 2);

y = (A/B - p).^2;

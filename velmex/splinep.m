function p = splinep(per,amp)
% splinep:  spline smoothing parameter for call to Spline Toolbox csaps
% p = splinep(per,amp);
% Last revised 2005-6-1
% 
% Spline smoothing parameter for call to Spline Toolbox csaps
%
%*** IN 
%
% per (1 X 1)r - period in years at which amplitude of freq response should be amp
% amp (1 x 1)r - desired amplitude of freq response at period per years
%
%*** OUT 
%
% p (1 x 1)r spline smoothing parameter
%
%*** REFERENCES 
% 
% I do not yet have a published reference for the equation below.  The equation is 
% presented courtesy of Jean-Luc Dupouy.
%
%*** NOTES
%
% The equation below is also described in the time series  notes and is given there as eq 8. 
%
% p =   1  /  (1 + (1 - u(f))*(cos(2 * pi * f)+2) / (12 * u(f) * (cos(2 * pi * f) - 1)^2))
%
% f = 1/per: target frequency (e.g., f=0.01,  corresponding to wavelength 100 yr)
% u(f) == amp: desired amplitude of frequency response at target frequency  f
%
% The spline parameter compute by the formula above is NOT identical to the spline parameter p
% described by Cook and Peters (1981).  Jean-Luc Dupouy pointed this out in an email to Ken 
% Peters (cc'd to me) on 2005-5-27. The spline parameter p by this function is appropriate if
% Matlab Spline Toolbox function csaps is to be used to generate the spline with the
% desired degree of smoothing. 
%
% This function has been checked to the extent that the correct amplitude damping is given
% when function csaps is used to generate the smoothed curve for a sythetic input time series
% that is a pure sinusoid.  Settings of   10<per<1000 and   0.05<amp<0.95 were used in the tests
% with satisfactory results.
% 

f=1.0/per; % period to frequency in cycles per unit of time

p = 1 / (    ( ( cos(2*pi*f) + 2 )*(1-amp) / (12 * amp* ( cos (2*pi*f) - 1 )^2 ) ) + 1    ); % Spline parameter




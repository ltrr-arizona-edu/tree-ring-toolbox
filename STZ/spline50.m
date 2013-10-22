function per=spline50(p)
% spline50(p): % find period with 50% frequency response for spline with parameter p
% per=spline50(p);
% Last revised: 9-9-99
%
% Called by crvfit.m if needed to get period where spline freq response is 50%
%
%*** INPUT
% 
% p (1 x 1)r spline parameter, after Cook and Peters, 1981
%
%*** OUTPUT
%
% per (1 x 1)r period (years) at which ampl of frequency response is 0.5 


freqstart=.01;  % starting frequency for iteration (period of 100 yr)

% Call subfuncton fun1 to iterate to 50% spline period
[f50,fval,exitflag]=fminsearch('splmin01',freqstart,[],p);
per = 1/f50;


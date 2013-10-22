function [g,w,p,eflag]=splpos(x,xcrit,t,tt,delta);
% splpos:  find spline such that smoothed series (growth curve) not negative or zero anywhere
% [g,w,p,eflag]=splpos(x,xcrit,t,tt,delta);
% Last revised 2-22-02
%
% Find spline such that smoothed series (growth curve) not negative or zero anywhere.  For time series
% of length N, start with an 2*N-yr cubic smoothing spline. If that smooth line is not everywhere above the
% critical level xcrit, fit a spline of wavelength 2*N/(1 + delta). If that curve still not above xcrit everywhere,
% try spline of wavelength 2*N/(1+2*delta), etc.  Thus splines of increasing frequency are tried until the smooth
% curve g is above xcrit everywhere.  
%
%*** INPUT
%
% x(mx x 1)r time series before smoothing
% xcrit (1 x 1)r  threshold above which the spline curve must fall
% t (? X 1)i   time vector for x
% tt (? x 1)i specified time vector for the output g (see notes)
% delta (1 x 1)r increment added to denominator ( see above in making spline more flexible) (see notes)
%
%*** OUTPUT
%
% g (? x 1)r  spline-smoothed version of x, at time values tt
% w (1 x 1)r  50% wavelength of the final spline (see notes)
% p (1 x 1)r  spline parameter for the selected spline with 50% wavelength w
% eflag (1 x 1)i error flag
%   ==0 no error
%   ==1 spline blows up: var(g) greater than var(x). See notes. 
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED 
%
% cfspl -- curve fit by spline
%
%*** TOOLBOXES NEEDED 
%
% Spline
%
%
%*** NOTES
%
% 50% wavelength.  Wavelength at which the amplitude of frequency response of the spline is 0.5
%
% Spline blows up.  Trials with csaps function have shown that the spline becomes oscillatory with high variance
% when the paramter gets very large.  This corresponds to a short wavelength, about 5 yr.
%
% delta:  I have used delta=0.2 with good results with lfcheck.m
%
% tt:  usually, this is the same as t

amp=0.5;
if size(x,2)~=1 | size(t,2)~=1 |   size(tt,2)~=1
    error('x, t and tt must be col vectors');
else;
    nx=length(x); % length of time series
end;
if length(t) ~=length(x);
    error('x and t must be same length');
end;
if delta <=0 | delta>=1;
    error('delta must be greater than zero and less than 1');
end;
if size(xcrit,1)~=1 | size(xcrit,2)~=1;
    error('xcrit must be scalar');
end;



% Fit increasingly flexible spline, beginning with 1*N, where N is length of series
kwh5=1;
ptry=2*nx; % initially try wavelength twice series length
denom=1; % denominator of factor setting the spline wavelength
while kwh5; % while the growth curve has a negative value
    ptry=ptry/denom;
    pptry=splinep(ptry,amp); % spline parameter
    g=(cfspl(pptry,tt,t,x))';
    L=tt>=t(1) & tt<=t(end);
    if sum(L)~=length(t);
        error('From splpos, line 99: tt does not include all of t');
    else;
        vratio=var(g(L))/var(x);
    end;
    if vratio>1; % Variance of spline greater than variance of original series (spline blew up)
        eflag=1;
        w=[];
        p=[];
        g=[];
        kwh5=0;
    else;
        if any(g<=xcrit);
            denom=denom+delta;
        else;
            % g acceptable
            w=ptry; % 50% wavelength
            p=pptry; % apline parameter
            kwh5=0;
            eflag=0;
        end;
    end;
end; % while kwh5; % while the growth curve has a negative value
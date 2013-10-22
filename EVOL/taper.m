function xt = taper(x,p)
%
% Split-cosine-bell tapering of a time series.
%
% D. Meko, 7-2-95
%
%***************** IN ARGS ******************************
%
% x (mx x 1)r  time series
% p (1 x 1)r   total proportion of series to be tapered
%		(1/2)p of each end will be tapered
%
%
%********************* OUT ARGS **************************
%
% xt (mx x 1)r taperered time series
%
%
%
%********************* NOTES **************************
%
% Source: Bloomfield  1976, . 116
%
% Time series x should have been converted to zero-mean (e.g., by
% detrending) before calling taper.m
%


% Size and check arguments
[mx,nx]=size(x);  % mx is length of time series (N in Bloomfield)
if mx<20 | nx ~=1,
	error('x must be col vector of minimum length 20')
end
if p<=0 | p>=1; % taper proportion is capital P in Bloomfield
	error('Taper proportion must satisfy 0<= p <=1')
end


M=(round(p*mx))/2;  % number of values on each end of series that will
		% be modified

I = (1:M)';

weight= 0.5 - 0.5*(cos(pi*(I-0.5)/M));



% Taper first part
x(I) = x(I) .* weight;

% Taper last part
Itail = ((mx-M+1):mx)';
x(Itail) = x(Itail) .* flipud(weight);


xt=x; % return argument, the tapered version of input x

function [N,m,ps,p]=binom1(x,y,k,q)

% binomial test of number of successes m in N trials.
% See Panofsky and Brier, p. 33.
%
%  D. Meko , 12-17-93
%
% 
% *********   INPUT  *******************
%
% x (nx x 1) one time series, nx observations
% y (ny x 1) another time series, ny=nx observations
% k (1 x 1) option control
% 
%		k=1 "below qth quantile":  elements of x and y below their
%		    respective qth quantiles are marked.  The values of y
%		    in the marked years of x are considered the "sample".
%		    A success is a sample value for which y is also below its
%		    qth quantile.
% q (1 x 1)  the quantile to be used as threshold. E.g., q=0.1 means
%		the 0.1 quantile -- or in other words the lowest decile.
%
%
%**********************  OUTPUT  ************************
%
% N (1 x 1) the number of "samples", defined as above as the 
%	number of x values below the threshold
% m (1 x 1) the number of successes.  A success is a sample value
%  for which y also is below its threshold
% p (1 x 1) the computed probability of observing m successes in 
%   N trials



if k==1;  % So far the only option -- the below quantiles option
	Lx = x <=quantile(x,q);
	Ly = y <=quantile(y,q);
	N=sum(Lx);  % number of samples
	ps = sum(Ly)/ length(y);  % prob that a value of y is below its
%			threshold
	m = sum(Lx & Ly);   % number of successes

else
end


	p = ((prod(N:-1:1))/((prod(m:-1:1)) * (prod((N-m):-1:1))) ...
		* (ps ^m) * (1-ps) ^(N-m));
	

	

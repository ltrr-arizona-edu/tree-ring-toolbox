function [Nsamp,nhit,p]=hypernew(x,y,k,q)
% hypernew:  hypergeometric test for coincidence of values outside threshold
% CALL: [Nsamp,nhit,p]=hypernew(x,y,k,q);
% 
%  Meko 4-8-99
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
% Nsamp(1 x 1) the number of "samples", defined as above as the 
%	number of x values below (if k==1) the threshold
% nhit (1 x 1) the number of successes.  A success is a sample value
%  for which y also is below its threshold
% p (1 x 1) the computed probability of observing "at least m" successes in 
%   Nsamp trials


if k==1;  % So far the only option -- the below quantiles option
	Lx = x <=quantile(x,q);
   Ly = y <=quantile(y,q);
   Nsamp=sum(Lx);  % number of samples
   if Nsamp~=sum(Ly);
      error('same quantile give different number of exceedances in x and y');
   end
   Npotent=Nsamp; % Number of potential successes (possible hits in y) is equal
   %  to the number of samples chosen by small (or large) values in x
   
   nhit = sum(Lx & Ly); % number of successes 
   
   % Call to hygecdf will return probability of having nhit or fewer successes. But
   % want the prob that number of successes is nhit or greater.  Thus will call for
   % prob of nhit-1 or fewer successes and take 1 minus that probability.  But must
   % watch for case of no hits, because then nhit minus 1 is negative.  
   if nhit==0;
      pcum=[];
      p=1; % probability of having zero or more hits 
   else
      % Compute cumulative prob of fewer than nhit successes
      pcum = hygecdf(nhit-1,length(x),Npotent,Nsamp);
      % Compute prob of nhit or more successes
      p=  1 - pcum;
   end;
   
   

else
end
	

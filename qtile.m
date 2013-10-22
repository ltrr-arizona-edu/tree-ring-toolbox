function [q,f]=qtile(x,f)
% qtile: f-th quantiles of data
% [q,f]=qtile(x,f);
% Last revised 8-24-00
%
% Quantiles of a vector of data
%
%*** INPUT
%
% x (mx x 1)r   data
% f (mf x 1)r <optional>  f-values at which quantiles are to be computed. If
%   omitted, quantiles are computed at f-values of (1-0.5)/n, (2-.5)/n,...
%   (n-0.5)/n, where n==mx is the length of the data.  See notes contraints on 
%    f. 
%
%*** OUTPUT
%
% f (mf x 1)r  the f-values at which the quantiles were computed. If f was included
%   as an input argument, f on output is same as f on input. Otherwise, f is computed
%   as described under INPUT.
% q (mq x 1)r  the f-th quantiles of the data; the quantiles corresponding to the 
%   output f-values.  mf==mq.
%
%*** REFERENCES 
%
% Cleaveland, W. S., 1993.  Visualizing Data.  Hobart Press, Summit, New Jersey, 360 p.
%   [p. 16-20] 
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Algorithm:  rank data values from smallest (i=1) to largest (i=n). If f  not included
% as input argument, compute quantiles at f-values (i-0.5)/n, where n is the data length.
% If f included, compute as if f not included and interpolate linearly for desired 
% quantiles.
%
% f, if included as an input argument, must be monotonic increasing no smaller than
% (1-0.5)/n and no greater than (n-0.5)/n, where n is the sample size

% Check data length
[mx,nx]=size(x);
n=mx; % n is now the sample size
i=(1:n)';
j=(i-0.5)/n;  

% Check for missing data in x
if any(isnan(x));
   error('x must contain no NaNs');
end;

% If f not specified as input arg, compute quantiles at each data point
if nargin==1;
   f=j;
end;

% Check  x and f dimensions
[mf,nf]=size(f);
if nx~=1 | nf~=1;
   error('x  and f must be column vectors');
end;

% Check that f monotonic increasing
if ~all(diff(f)>0)
   error('f must be monotonic increasing');
end;


% Cannot compute quantiles outside the data range
fmin=(1-0.5)/n; % minimum computable quantile
fmax=(n-0.5)/n; % maximum ...
if any(f<fmin) | any(f>fmax);
   error('data length insufficient to allow computation of minimum or maximum quantile');
end;

s=sort(x);  % Sort data;  s contains the j-th quantiles
q=interp1(j,s,f); % Interpolate f-th quantiles from j-th quantiles

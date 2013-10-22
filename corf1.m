function R=corf1(x,y,M)

% Correlation matrix corresp to covf.m crossproducts mtx
% Works for two time series only
%
% D. Meko  7-13-93
% 
%**********  INPUT ARGS  *************************
% 
% x (mx x 1) first time series, usually the master series
% y (my x 1) second times eries, usually the sample (undated)
% M number of + and  - lags to compute xcorrs
%
%
%***********  OUTPUT ARGS  **************************
%
% R matrix of correlation coefs, as follows
%
%  rows:	1 - lag k autocorrelations of x
%		2 - cross corrs of y lagged 0,-1,-2,..-M 
%			from x
%		3 - cross corrs of y lagged 0, 1, 2,... M
%			from x
%		4 - lag k autocorrels of y
%
%  cols:  col 1 is lag-0.  Col 2 is lag 1.  ...
%		col M is lag +- lag (M-1)
%
%***********   NOTES  *********************************
%
%  x and y must be column vectors



[mx,nx]=size(x);
[my,ny]=size(y);

if (nx~=1) | (ny~=1), error('x and y must be col vectors'),end
if (mx ~= my), error ('x and y must be same length'), end 
  
if (mx/M < 4), disp ('You set number of lags to > 1/4 n of obs'),end

z=[x y];
zs=std(z) * sqrt((mx-1)/mx);
zs=zs(ones(mx,1),:);

z1=dtrend(z);  % subtract col means
z2=z1./zs;  % convert to z-scores

R=covf(z2,M+1);  % crossproducts matrix of z-scores = correls

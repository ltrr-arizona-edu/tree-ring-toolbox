function [N,m,ps,n95,n99,p]=hygeom1(x,y,kq,ks,klq)
%
% USAGE : [N,m,ps,n95,n99,p]=hygeom1(x,y,kq,ks,klq)
%   Hypergeometric test of number of successes m in N trials.
%
%  S. Anwar	8/4/94.
%
%  
% INPUTS 
%-------
% x (mx x nx) 	Time series, nx observations
% y (my x ny) 	Another time series, ny=nx observations
% kq (1 x 1)	If kq=1, perform lower 1-sided test
%		If kq=2, perform upper 1-sided test
%		If kq=3, perform 2-sided test
% ks (1 x 1)	If ks=1, typical RW mode
%		If ks=2, typical skeleton mode
%		If ks=3, mixed mode
%
%
% OUTPUTS  
%--------
% N (nx x 1) 	The number of "samples", defined as above as the 
%		number of x values below the single-sided threshold 
%		or outside the double-sided threshold
% m (ny x 1) 	The number of successes.  A success is a sample value
%  		for which y also is below its single-sided threshold
%		or outside the double-sided threshold
% ps (1 x ny) 	The number of y's outside qu and/or ql percentile
% n95 (1 x ny)	95 % significance level
% n99 (1 x ny) 	99 % significance level
% p (1 x ny) 	Computed probability of observing m successes in 
%   		N trials
%
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% QUANTILE.M	Picks the data values which are outside the threshold
% HYGHLP.M	Computes the Hypergeometric statistics
% SLVMNU.M	A modified menu function
%________________________________________________________________
 
% ql (1 x 1)	Lower threshold quantile. E.g., ql=0.1 means
%		the 0.1 quantile -- or in other words the lowest decile.
% qu (1 x 1)	Upper threshold quantile

[my,ny] = size(y);
if ks==1,	% RW or IND mode - Both x and y contain RW or IND data
  if kq==1,
    mnul=['0.1';'0.2';'0.3';'0.4';'0.5';'0.6';'0.7';'0.8';'0.9'];
    ql=slvmnu('Enter LOWER ql',mnul);
    ql=0.1*ql;
    qx = quantile(x(:,1),ql); 
    qxl = qx*ones(1,ny);	% qxl (1 x ny)
    qyl = quantile(y,ql);   	% qyl (1 x ny)
    qxml = qxl(ones(1,my),:);
    qyml = qyl(ones(1,my),:);
    Lxl = x <= qxml;
    Lyl = y <= qyml;
    [m,N,ps,n95,n99] = hyghlp(my,Lxl,Lyl,klq);
  elseif kq==2, 
    mnul=['0.1';'0.2';'0.3';'0.4';'0.5';'0.6';'0.7';'0.8';'0.9'];
    qu=slvmnu('Enter UPPER qu',mnul);
    qu=0.1*qu;
    qx = quantile(x(:,1),qu); 
    qxu = qx*ones(1,ny);	% qxl (1 x ny)
    qyu = quantile(y,qu); 	% qyu (1 x ny)
    qxmu = qxu(ones(1,my),:);
    qymu = qyu(ones(1,my),:);
    Lxu = x >= qxmu;
    Lyu = y >= qymu;
    [m,N,ps,n95,n99] = hyghlp(my,Lxu,Lyu,klq);
  else
    error('kq must be either 1 or 2');
  end
elseif ks==2,		% Skeleton mode - Both x & y skeleton
  qm=zeros(my,ny);
  if kq==1,  
    Lxl = x < qm;	% Ones if x below zero
    Lyl = y < qm;	% Ones if y below zero
    [m,N,ps,n95,n99] = hyghlp(my,Lxl,Lyl,klq);
  elseif kq==2,
    Lxu = x > qm;	% Ones if x above zero
    Lyu = y > qm;	% Ones if y above zero
    [m,N,ps,n95,n99] = hyghlp(my,Lxu,Lyu,klq);
  else
    error('kq must be either 1 or 2');
  end
elseif ks==3,	 	% Mixed mode - x skeleton, y ringwidth or index
  qm=zeros(my,ny);
  if kq==1,
    Lxl = x < qm;		% Ones if x below zero
    ql = sum(Lxl)/my;		% Proportion of x-values below zero
    qyl = quantile(y,ql(1));	% Equivalent lower quantile of y
    qyml = qyl(ones(1,my),:);
    Lyl = y <= qyml;	% Ones of y below threshold
    [m,N,ps,n95,n99] = hyghlp(my,Lxl,Lyl,klq);
  elseif kq==2,
    Lxu = x > qm;
    qu = ones(1,ny)-sum(Lxu)/my;	% Sum of x-values above zero
    qyu = quantile(y,qu(1)); 		% Equivalent upper quantile for y
    qymu = qyu(ones(1,my),:);   
    Lyu = y >= qymu;
    [m,N,ps,n95,n99] = hyghlp(my,Lxu,Lyu,klq);
  else
    error('kq must be either 1 or 2 ');
  end
else
    error('ks must be either 1, 2, or 3');    
end

if nargout==6,
   
  % Calculate the probability using hypergeometric distribution

  [mxm,mxi] = max(m);
  p = hygecdf(mxm,my,ps(mxi)*my,N(1));

elseif nargout==5;
else
  error('Invalid number of outputs');
end

% End of file

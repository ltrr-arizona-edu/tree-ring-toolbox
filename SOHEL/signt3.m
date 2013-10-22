function [na,ns,nkey]=signt3(x,y,k)
%
% USAGE : [na,ns,nkey]=signt3(x,y,k)
%   Given two time series matrices of the same size : say X and Y
%   Each column of X is a time series.
%   Each column of Y is a time series.
%   This function does a series of sign tests; column 1 of X with 
%   column 1 of Y, 2 with 2 and so on.
%   Test two time series for like sign departures from respective means,
%   or for like sign of first differences.
%
%
% INPUTS
%-------
% x (my x ny)	A time series matrix
% y (my x ny)	Another time series matrix
% k (1 x 1)	k = 1 : test for like sign of departures from respective means
%		k = 2 : test for like sign of first difference
%
%
% OUTPUTS 
%--------
% na (? x 1)	na(1) is the number of sign agreements of X(:,1) with Y(:,1) 
% ns (? x 1)	ns(1) Total number of cases for tests for cols of X vs Y. 
%               Might be less than the number of rows in X and might vary
%               from column to column.
% nkey(:,1)	If n is significant at 80% level
% nkey(:,2)	If n is significant at 90% level
% nkey(:,3)	If n is significant at 95% level.
% nkey(:,4)	If n is significant at 99% level.
%
%
% USER WRITTEN FUNCTIONS NEEDED 
%------------------------------
% SGNTBL2.TAB	A table of 95% and 99% significance levels. Must be 
%               available in matlab search path.
%_____________________________________________________________________
  
  [m1,n1] = size(x);
  [m2,n2] = size(y);
  L = [n1~=n2 m1~=m2];
  if any(L);
    error('Input error: x and y must be of same size');
  end

  % Subtract means from cols of x, y

  xmean = mean(x);
  ymean = mean(y);
  xmean = xmean(ones(m1,1),:);
  ymean = ymean(ones(m2,1),:);
  x     = x - xmean;
  y     = y - ymean;

  if k==1;                      % Test departures from means
    xs = sign(x);               
    ys = sign(y);
  elseif k==2;                  % Test first-differences
    xs = sign(diff(x));
    ys = sign(diff(y));
  else
    error('signt3.m: k must be either 1 or 2');
  end

  ns    = m1 * ones(n1,1);      % Vector of unadjusted sample sizes
  L2    = xs==0 | ys==0;        % Logical pointer to "ties" in sign comparison
  ns    = ns - sum(L2)';        % Adjust for ties

  if ns < 8 ;
    disp(ns);
    error('fewer than 8 cases');       % See Fritts(1976)
  end;

  L     = xs==ys & xs ~= 0 & ys ~= 0;  % Match in sign
  na    = sum(L)';              % row vector : number of agreements in each comparison
  nd    = ns-na;                % number of disagreements
  n     = min([ na'; nd']);     % The test statistic, row vector.
  n     = n';                   

  %     Look up table value to test for 95%  and 99% significance, or
  %     ns >= 50, compute the relevant values. Ignore ties in dealing to 
  %     use lookup table or the equation.

  z=[1.28 1.64 1.96 2.58];
  if max(ns) >= 50;                  % Compute significant levels.
      nkey = [(ns-1-z(1)*sqrt(ns))/2, (ns-1-z(2)*sqrt(ns))/2,...
              (ns-1-z(3)*sqrt(ns))/2, (ns-1-z(4)*sqrt(ns))/2];
  else
    if ~exist('sgntbl2');
      load sgntbl2.tab;
    end
    nkey = table1(sgntbl2,ns);
  end
  nkey=max((ns(:,ones(1,4))-nkey),nkey);
end             

% End of the function

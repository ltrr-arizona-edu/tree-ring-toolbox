function [T,I0,T0,dstari0,sig,T05]=potter(x,y,tablepot)

% Test for detecting shift in mean of a precipitation series

% Reference: "Illustration of a new test for detecting a shift in mean
% in precipitation series", by Kenneth W. Potter, Monthly Weather
% Review 109, 2040:2045.

%************   INPUT ARGUMENTS   *******************************

% x (m1 x 1) precipitation time series thought to be "good"--
%		homogeneous.
% y (m1 x 1) precipitation series to be tested for homogeneity.
% tablepot (8 x 5) table of critical values (Potter, 1981, p. 2041)

%***********   OUTPUT ARGUMENTS   *******************************

% T    (m1-1 x 1) corresp to  T   in Potter, p. 2041.  m1 is the
% 		number of observations in y.  
% I0   (1 x 1) row index pointing to the year of the largest value of T
% T0   (1 x 1) the test statistic -- largest value of T
% dstari0   the amount  of change indicated in the mean of series y.
%		Meaningful only if T0 turns out to be significant.
%		To adjust y,  dstari0 should be added to y beginning with the
%               first year after that pointed to by I0.
% sig  (1 x 1) significance level of T0. Can take 5 possible values:
%     .01 -- most significant
%     .05
%     .10
%     .25
%     []   -- sig level greater than .25 (least significant class)
% T05 (1 x 1) value of T corresponding to sig level of .05.  T05 is
%		useful for plotting along with a time series of T against year
% 		to show when the inhomogeneities arise.


%***************   NOTES   **************************************

% Assumptions:
%
% 1. {x y} an independent sequence of two-dimensional vectors
% 			 each distributed bivariate normal
% 2. sequence {x y} is stationary, with the exception of a possible
%    shift in mean of y
%
%  Potter discusses the validity of these assumptions for precipitation 
%  data, and the consequences of violating the assumptions.
% 
%  The is no reason that the method cannot be used for any other type of 
%  time series data that satisfies the assumptions (perhaps monthly mean
%  temperature).

%  Maximum allowable sample size = 100, because table
%  of critical values in Potter (1981) has not been developed for 
%  larger sample sizes.

%  Be sure to load tablepot into the workspace before calling this
%  function.

%**************   BEGIN CODE   *************************************


% Check number of arguments
	if nargin ~=3, error('WRONG NUMBER OF INPUT ARGS'), end;
	if nargout ~=6, error('WRONG NUMBER OF OUTPUT ARGS'), end;


if ~exist('tablepot')
	load tablepot
end
 


m1=length(x);
m2=length(y);
if m1~=m2
	error('x and y are of different length')
end


% To use significance test, must standardize x and y.  Let u and 
% v be the corresponding standardized variables, with mean 0 and
% stdev 1.

xm=mean(x);
ym=mean(y);
xstd=std(x);
ystd=std(y);

u=(x-xm(ones(m1,1),:)) ./ xstd;
v=(y-ym(ones(m1,1),:)) ./ ystd;


%******   Begin coding Potter's equations   ******

I=(1:m1)';
X=cumsum(u) ./ I;
Y=cumsum(v) ./ I;

Xbar=X(m1);   Ybar=Y(m1);


s1x=u-Xbar;    s1y=v-Ybar;
Sx = sum(s1x .* s1x);  Sy = sum(s1y .* s1y);
Sxy = sum(s1x .* s1y);

N=m1(ones(m1,1),:);


X1=X-Xbar;    X2=Xbar-X;
Y2=Ybar-Y;

Sx(m1)=[];  Sy(m1)=[];  Sxy(m1)=[];  N(m1)=[];  I(m1)=[];
Y2(m1)=[];  X2(m1)=[];  I(m1)=[];  X1(m1)=[];

NmI=N-I;


F = Sx- (X1 .* X1 .* N .* I) ./ NmI;
D = (Sx .* Y2 - Sxy .* X2) * m1 ./ (NmI .* F);
T = (I .* NmI .* D .* D .* F) ./ (Sx .* Sy - Sxy .* Sxy);


% Find the maximum value of T, store it and the corresponding
% amount of change dstari0.  Factor ystd comes into play because
% want dstari0 in terms of original units of y, not of standardized
% units.


[T0,I0]=max(T);  % I0 holds the row index of max value of T.

dstari0 = D(I0) * ystd;


% Find the significance level of T0.  Only the discrete levels
% .01 .05, .10, and .25 are considered.  If sig level > .25, 
% sig is set to [].

t=[table2(tablepot,m1,.25)  table2(tablepot,m1,.10) ...
table2(tablepot,m1,.05)    table2(tablepot,m1,.01)]';

if T0< t(1)
	sig=[];
elseif T0>= t(1) & T0 <t(2)
	sig=.25;
elseif T0 >= t(2) & T0 < t(3)
	sig=.10;
elseif T0 >= t(3) & T0 < t(4)
	sig=.05
elseif T0 > t(4)
	sig=.01
end

T05 = t(3);  % Critical value corresponding to sample size m1 and
%		signif level .05


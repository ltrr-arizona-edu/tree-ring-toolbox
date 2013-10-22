function [na,ns,n95,n99]=signtest(x,y,k1)
% Sign test for agreement between two time series (p. 330, Fritts 1976)
% Test two time series for like sign departures from respective means, 
% or for like sign of first differences.
%
% Last revised 10-23-02 
%
%******   INPUT ARGS
%
% x (m1 x 1)  a time series
% y (m2 x 1) another time series
% k1 (1 x 1)  = 1   test for like sign of departures from respective means
%             = 2   test for like sign of first difference


%*****   OUTPUT ARGS
%
% na (1 x 1)  number of agreements in sign
% ns (1 x 1)  total number of cases for test.  Might be less than number
%		of years in time series.
% n95 (1 x 1) T or F:  T= n is significant at 95 % level
%			         F= n is not significant at 95 % level
% n99 (1 x 1) T or F:  like n95, except for 99% level


%*******  OTHER USER FUNCTIONS NEEDED
%
%  signtbl.mat  - table of 95% and 99% significance levels;  must be
%		available in matlab path


%*****   Size arrays; test that x and y are row vectors;  test that x,y same
% length;   subtract sample means from x,y;  Initialize significance-test
% results as "false"

[m1,n1]=size(x);  [m2,n2]=size(y);
L1 = [n1~=1  n2~=1  m1~=m2];
if any(L1);
	error('x and y must be CVs of same length!');
end

x=x-mean(x);
y=y-mean(y);

n95='N';  n99='N';


%*********  Optionally test departures or change
% Convert each series to either 1, -1 or zero

if k1==1;  % departures
	xs=sign(x);
	ys=sign(y);
elseif k1==2
	xs=sign(diff(x));
	ys=sign(diff(y));
else
	error('  k1 must = 1 or 2 !');
end


%******  Compute number of matches of 1's and -1's
% First eliminate cases where either xs or xy is 0.

L2= xs==0  | ys==0;
xs(L2)=[];   ys(L2)=[];
ns=length(xs);  % number of valid cases for test
if ns < 8;  % Fewer than 8 cases, see Fritts (1976)
	disp('Warning: fewer than 8 cases ');
end

L3=  xs==ys;  % match in sign
na=sum(L3);  % number of cases with agreement in sign
nd = ns-na;  % number of disagreements
n = min([na nd]);  %test statistic

%*****  Look up table value to test for 95% and 99% significance, or
%	if ns >= 50, compute the relevant values.

if ns >= 50;  % Compute significant levels
	nkey = [(ns-1-1.96*sqrt(ns))/2   (ns-1-2.8*sqrt(ns))/2];
else;  % use lookup table 
	if ~exist('signtbl');
		load signtbl;
	end
	%nkey=table1(signtbl,ns); % OBSOLETER
    nkey=interp1(signtbl(:,1),signtbl(:,2:size(signtbl,2)),ns); %  REVISED 10-23-02
end

if n <= nkey(1), n95='Y'; end;  % sig at 95%
if n <= nkey(2), n99='Y'; end;  % sig at 99%

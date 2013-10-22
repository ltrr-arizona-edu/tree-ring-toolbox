function y2=freqest1(x1,y1,x2)
% freqest1.m  missing value by cum freq distrib estimation
%
%
% x1 -- cv of sample of values from predictor station to be used to build the
%	freq distribution for that series
% y1 -- cv of sample of values from key station, for similar use
% x2 -- cv of values from pred station for those observations with missing data
%	at the key station
%
%
% y2 -- cv of estimated values for the key station for the same years represented by
%	x2
%
%
%********   USER - WRITTEN FUNCTIONS NEEDED -- NONE
%
%
% D. Meko 10-30-93
%
% Low-level function for use in estimating missing daily precipitation values for a
% group of stations in NPS SODE project.  General Idea:
% 
% Need to estimate data for a "key" station 
% Have daily ppt for the days in question at a predictor station
% If zero ppt at pred station, set ppt at key station to zero
% Otherwise, using the cum freq distributions based on only
%   non-zero values, interpolate the missing value as the value with the same 
%   empirical exceedance prob as the value for the same date at the pred station
% This function intended to be called separately for each of 366 Julian days
% Handled in calling functions is culling out of appropriate rows from matrices to
%	pass as input arguments.  In particular, the frequency distributions are to be built
%	from 3-day groupings.  For example, computation of estimates for Julian day 255 will
%	receive x1 and y1 samples based on all non-missing daily ppt values for Julian days 254,255,256
%	This allows larger sample size for freq dist computation, and smooths over period not longer than
%	typical weather systems -- so freq dist should be representative of middle day.



%*********  GLOBALS -- NONE


mx1=length(x1);  % sample size for cum freq dist of pred series
my1=length(y1);  % sample size for cum freq dist of key series
mx2=length(x2);  % number of daily values to be estimated


% Cum cum prob dist for samples x1 and y1.  Make lookup table from x1 with the
% prob as col 1 and values as col 2.  Make table for y1 with the values
% as col 1 and the probs as col 2.
sx1=sort(x1);
% Need small increment to ensure that sorted x1 is monotonic increasing
% Find smallest x1; divide by 100 to get 1% of it;  Make increment
% go from minus to plus 1 percent of smallest x1
inc=(linspace(-sx1(1)/100,sx1(1)/100,mx1))';

sx1=[sx1+inc  ((1:mx1)/(mx1+1))'];  % exceedance prob in col 1
sy1 = sort(y1);
sy1 = [((1:my1)/(my1+1))'  sy1];  % exceed prob in col 2

% Lookup exceed prob assoc with each value in x2
p = table1(sx1,x2);


% Check that no values in p fall outside the upper and lower cum prob values
% in the cum prob table for the key series.  If any such values occur, replace
% with the corresponding lowest or highest observed values in y1
min1=min(sy1(:,1)); %minimum tabulated exceed prob in table for key series
max1=max(sy1(:,1)); % max tabulated excced prob
ilow=find(p<min1); % index to rows of p below lowest accept prob
ihi =find(p>max1); % index to rows of p below highest accept prob
p(ilow)=min1(ones(length(ilow),1),:);
p(ihi)=max1 (ones(length(ihi),1),:);


% Interpolate missing values for key series
y2 = table1(sy1,p);

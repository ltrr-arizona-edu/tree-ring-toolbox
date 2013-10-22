function [Y,Z,B,R]=monfill(X,W)

% Fill in missing monthly precipitation or temperature values by
% median ratio or median departure method


% D. Meko, begun 2-24-94

% X (mX x nX) matrix of monthly ppt for nX-2 stations and 
%  tota months.
%	Col 1 is the calendar year.  Col 2 is the month (Jan=1).  Remaining
%	cols are the ppt or T values for each station. 
% W (mW x 5) pointer to cols of X holding predictor stations.
%	Each row stands for a key station, mW total stations.
%	Col 1 gives the sequence number in X of the most preferred predictor 
%	station.  W is built before running this function by common sense by 
%	setting priorities of stations based on proximity and similar climate
%	regimes.  After running program, spearman r and other output gives 
%	feedback on how reasonable your priorities are.
%
%
% Y (mY x nY)  the filled-in equivalent of X.  Should contain no missing
%	values.
% Z (12 x nS) for each monthy, the proportion of days with zero
%	precipitation.  Based on data before estimation.  nS=nW is the number
%	of stations.
% B (366 x nS) p-value that observed number of successes in binomial test
%	could occur by chance.  Based on comparison of data for each key
%	with its highest-priority predictor station only.
% R (366 x 5xnS) spearman R between key stations and each of its 5 pred
%	stations.  Each row for a Julian day.  Cols 1-nS for the top-priority
%	predictor station.  Cols nS+1 to 2nS for the second-priority, etc.
%	Spearman coef computed only on those observations with nonzero ppt at
%	both stations.
%
%
%*********    METHOD *************
%
% Given a key station with some missing daily ppt values, and a set of
% other stations with daily ppt data.  Using climatological reasoning, set
% a priority for the other stations as predictor stations
% for estimating missing data at the key station.  The priority can be
% based on distance between stations, similarity of climatological
% or some other criterion.
%
% Take the highest-priority predictor station.  Assume that a missing 
% daily ppt value at the key station will have the same empirical 
% exceedance probability (disregarding zero-ppt values) as the observed
% daily ppt value at the predictor station. For a given Julian day, the
% samples for computing the frequency distributions are values for that
% same Julian day and the two flanking Julian days.  Say Julian day 255
% in 1952 is misssing from the key station.  The empirical prob 
% distribution of the pooled non-zero values of Julian-days 254,255,256 
% at the predictor station is computed.  The same for the prob distrib
% at the key station.  The distribution for the pred station is then used
% to interpolate the exceedance prob of the observed Julian-day 255 ppt
% in 1952 at the pred station.  The ppt value with the same exceed prob
% in the distribution for the key station is used as the estimate for
% the missing ppt.
%
% One predictor station might not have daily data for all the missing ppt
% at the key station.  Additional missing values are estimated as above,
% but using the second, third, ... priority predictor station.
%
% Because matrices can be huge for, say, 100 years of 366 days of ppt at
% 20 stations, this function is likely to drag on forever swapping out
% to disk for some datasets, especially if RAM is low.  So sometimes it
% might be best to manually carry out the functions of this highest-level
% function separately by invoking subroutines from the workspace. 
% The idea would be to deal with each pair of key/predictor stations 
% separately to get a new "estimated" column for the key station. 
%
% Nesting of functions:

%  monorg.m -- before running monfill.m, use monorg.m to reorganize
%		monthly matrix
%	monfill.m -- top level function; has loop over key station, and
%		   inner loop over predictor stations; counts number of 
%		   missing values, zero values, etc for each Julian day at
%		   each station
%		








% 		pairchk.m -- counts number of missing values for each Julian
%		   day and station.  Calc pctg of non-missing values that
%		   is zero for each Julian-day and station.
%		  
%		dayest.m  -- loops through 366 julian days; handles bookkeeping
%			tasks with subscripts so that vector of estimates easy
%			to overwrite modified mother array; handles special cases
%			of no estimates needed or possible, and of zero-value
%			estimates
%
%			binom.m -- binomial computations; what is pvalue that
%			  observed number of successes (rainy day at key stn
%			  for rainy day at pred stn) occurs by chance
%
%			spear.m -- spearman correlation coefficient.  
%			  Considering nonzero ppt values only, does higher ppt
%			  at pred station go with higher ppt at key stn?
%
%			freqest1.m -- computes freq distribs and estimate 
%			  missing daily ppt values.  The "meat" function of the
%			  whole procedure.





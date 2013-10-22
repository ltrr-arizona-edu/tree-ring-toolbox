function [P,Z]=pairchk2(X);
% Proportion of daily ppt values missing or zero for various Julian days
% and stations -- fully vectorized version
%
% X (mX x nX) matrix of daily ppt at several stations.  Col 1 assumed to be
%	year,  col 2 the Julian day, remaining cols the ppt for several stns
%	X must cover 366 days for each of however many years of data there are.
%
% M (366 x nX-2) number of missing values for each Julian day at each stn
% Z (366 x nX-2) proportion of non-missing values that are zero
%
%
%******* D. Meko   11-1-93
%
%
%*****  CONTEXT  ************
%
% One of many functions under top-level function dayfill.m, used to estimate
% missing daily ppt in a matrix of ppt data
%
%
%******** NOTES
%
% Missing values assumed to be coded as less than zero ppt
%
%
%*********  GLOBALS -- none
%
%

[mX,nX]=size(X);
nyr=mX/366;  % number of years of data
nS= nX-2;  % number of stations

U=zeros(nyr*366,1);
L=zeros(nyr,366*nS);
P=zeros(366,nS);
Z=zeros(366,nS);

% Say have 100 yr of data and 15 stations.  Then mX=36,600 and
% nX = 17.  These example data sizes will be used in the discussion
% below.

% The general idea in vectorizing is to avoid loops.  I do this by building
% massive arrays of subscripts, and pulling out desired data using the
% subscript arrays.

% First step is to build matrix to pull out the Julian days 1-366.
% The row indices of day 1 in X are 1, 367, 367+366, etc.  The indices for
% day 2 are 2, 368, 368+366, etc.  And so forth through day 366.  Build
% matrix whose first col has subscrpts for day 1, second col for day 2, and
% so forth. 
x1=(1:366:366*nyr)';  % cv of subscripts for Julian day 1
X1=x1(:,ones(366,1));  % dupe x1 to 366 identical cols
	% X1 is now of size nyr x 366

% Need to add numbers to X1 so that col 2 begins with 2, col 3 with 3,
% and so on.  Do this by building an nyr x 366 matrix of increments
% First make a row vector of 0 1 2 ...365.  Then dupe the row nyr times.
e1=0:365; % row vector
E1=e1(ones(nyr,1),:);  % e1 duped into nyr rows

% Add the increment array to the subscript array.
I1=X1+E1;  % I1 is (nyr x 366) subscript array to rows of X.

% String columns of I1 into a single column of row subscripts.
U(:)=I1;

% Need a col-subscript vector pointing to cols of X holding ppt values for
% each series.
V=[3:nX];

% Pull data for the selected rows and columns out of X
Y = X(U,V);  % For example, Y is (36,000 x 15)


%************  Tabulate Number of non-missing values for each Julian Day
%			at each station
%
LP=Y>=0;  % logical array with 1 for 'present' (non-missing) data,
%  zero for missing.  LP same size as Y.
	  
L(:)=LP;  % Re-shape LP columnwise to put nyr values corresp to Julian 
% day 1, station 1 into col 1 of L;  for day 2 into col 2, .... for
% day 366 into col 366.  Likewise, next 366 cols hold similar data for
% station 2, and so forth through all stations.  For example, 
% L is (100 yr x 366*15).

% Sum '1s' in a col to get number of non-missing values for that Julian
% day and station.
S1=sum(L);  % S1 is size (1 x 366*15)

% Re-shape S1 into P such that (1) columns refer to stations, and (2)
% rows refer to Julian days 1-366.  Each entry in P is the number of
% 'present' (non-missing) daily values for a given Julian day and station.
P(:)=S1';


%***********	Considering only the non-missing data, tabulate the
%			decimal fraction of data with zero ppt.  In other words
%			the fraction of the days with no rain for each Julian
%			day and station. 
%
% Logic in pulling out data and re-shaping of arrays is identical to that
% above, and is not repeated in comments.  
%
LZ=Y==0;
L(:)=LZ;
S1=sum(L);
Z(:)=S1';

% Convert Z (number of zero days) into a decimal fraction.  Round to
% nearest 3rd decimal digit (e.g., 0.925 means 92.5 percent of values
% are zero for a given Julian day and station).
Z=round((Z ./ P) * 1000) /1000;






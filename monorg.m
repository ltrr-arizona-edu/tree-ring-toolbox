function Y=monorg(X,P,c,xmiss)

% Organize monthly ppt or temp before running monfill.m
% 
% D. Meko, 2-24-94
%
%
%****************** INPUT ARGS  ****************
%
% X (mX x 13) sequential sets of "decks" of monthly ppt or temp
%		Col 1 is the year;  cols 2-13 are data for Jan-Dec
%
% P (mP x 2) starting and ending years for each monthly deck
% c (1 x 2) specified start and end year of output matrix Y
% xmiss (1 x 1) missing value numeric code 
%
%*****************  OUTPUT ARGS  ************************
%
% Y (mY x mP+2) reshaped version of X.  Col 1 holds year;  col
%   2 holds month;  cols 3:mP+2 hold the corresponding data
%	 for the mP stations.
%
%	 row  1 holds data for jan of year c(1)
%   row  2            for feb of year c(1) etc
%
%   row 13 hold data for jan of year c(1)+1 etc
%
%
%******************  USER-WRITTEN FUNCTIONS NEEDED  -- NONE
%
%
%*******************  NOTES *****************************
%
% You must make sure that all data in X are valid data or a valid
% numeric missing-value code, and that there are no completely
% missing rows (years) of data in any input deck.
%
% Fortran program 999fill.f comes in handy in the above


%***********  SIZE AND ALLOCATE

[mX,nX]=size(X);  
if nX~=13, error('X does not have 13 columns'), end;

[mc,nc]=size(c);
if ~all([mc nc]==[1 2]), error('c must be 1 x 2'), end;

[mP,nP]=size(P);
if nP ~=2,  error('P must have 2 columns'), end;



% Set a few index arrays and vectors into X
n=P(:,2)-P(:,1)+1;  % number of years in each station input record
b=sumsum(n);  % each station's record ends in this row of X
a = b-n+1;   %  eacg station's record starts in this row of X

% Set a few index vectors into Y
nc = c(2)-c(1) + 1;  % number of years in Y
mY= nc*12;  % number of rows in Y
nY = mP+2;  % number of cols in Y
Y=xmiss(ones(mY,1),ones(nY,1));   % intialize Y with missing values
e = (P(:,1)-c(1))*12 + 1;   % each station record starts in this row of Y
f = e + (n*12) -1;  % each station rec ends in this row of Y


% Loop by station to fill-er-up

for i = 1:mP;
	Xrows=(a(i):b(i));  % row range in X to get data from
	Yrows= e(i):f(i);  % row range in Y to stick data into
	X1=X(Xrows,2:13);  % pull out the block from X, minus the year col
	Y(Yrows,i+2)=reshape(X1',n(i)*12,1);  % string out as col, and subst
end

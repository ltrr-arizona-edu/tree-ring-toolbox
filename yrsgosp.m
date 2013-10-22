function P=yrsgosp(X,yr1);
% P=yrsgosp(X)
%
% Years of starting and ending non-NaN time series segments
%
% D Meko 12-23-95
%
%********************* IN ARGS **************************
% 
% X (mX x nX)r time-series matrix made up of continuous time series
%		with leading and/or ending NaN segments. Each col is a 
%		time series. No "years" column is included.Also works on a
%		vector time series, but must be column vector
% yr1 (1 x 1)i first year (row) or X
%
%********************** OUT ARGS ************************
%
% P (nX x 2)i starting and ending years of the period of valid
%		(not NaN) data for each column of X
%
%********************** UW FUNCTIONS NEEDED -- NONE **********
%
%******************** NOTES *********************************
%
% STRATEGY: Split X into three column-groups: leading NaNs, 
% trailing NaNs, and flanking NaNs.  Compute start and end years
% separately for the groups and re-combine years.
%
% Broken time series are not allowed.  In other words, NaNs can
% not be imbedded within the valid data of a column.  yrsgosp.m 
% checks for this problem.



% Size and allocate
[mX,nX]=size(X);

% Bombs on row vector
if mX==1, error('X must be matrix or column vector'), end;

a=NaN;
P=a(ones(nX,1),ones(2,1)) ; % will hold start, end years

% Compute end year of X
yr2 = yr1+mX-1;

% Identify columns
L=~isnan(X);
D=diff(L);
m=1-mX; % number of rows in D

% Check for broken time series
L1 = D==1 | D==-1;
if any (sum(L1)>2),
	error ('Imbedded NaNs within valid data a no-no')
end




%************** COLUMNS OF X WITH LEADING ZEROS *******
% Start year sometime after yr1, end year same as last year of X

c=any(D==1) & ~any(D==-1);
n = sum(c); % how many columns in this category
if n>0;
	i=find(c);  % column number(s) in X and P
	P1 = a(ones(n,1),ones(2,1));

	e1 = yr2(ones(n,1),:); % cv of ending years
	D1=D(:,i);
	s1= find(D1);
	seq1=(1:n)';
	P1(:,1)= s1+(m*(seq1-1)) + (yr1) ;  % start years
	P1(:,2) = e1;

	P(i,:) = P1;
end

%************** COLUMNS OF X WITH TRAILING ZEROS *******
% Start year in yr1, end year before yr2

c=any(D==-1) & ~any(D==1);
n = sum(c); % how many columns in this category
if n>0;
	i=find(c);  % column number(s) in X and P
	P1 = a(ones(n,1),ones(2,1));
	b1 = yr1(ones(n,1),:); % cv of beginning years
	P1(:,1)=b1;

	D1=D(:,i);
	s1= find(D1);
	seq1=(1:n)';
	P1(:,2)= s1+(m*(seq1-1)) + (yr1) -1;  % ending years
	P(i,:) = P1;
end


%************** COLUMNS OF X WITH LEADING AND TRAILING ZEROS *******
% Start year after yr1, end year before yr2

c=any(D==-1) & any(D==1);
n = sum(c); % how many columns in this category
if n>0;
	seq1=(1:n)';
	i=find(c);  % column number(s) in X and P
	P3 = a(ones(n,1),ones(2,1));
	
	DD=D(:,i);

	% Leading zeros first
	s1=find(DD==1);
	P3(:,1)= s1+(m*(seq1-1)) + (yr1) ;  % start years
	
	% Trailing zeros
	s1=find(DD==-1);
	P3(:,2)= s1+(m*(seq1-1)) + (yr1) -1;  % ending years
	P(i,:) = P3;
end

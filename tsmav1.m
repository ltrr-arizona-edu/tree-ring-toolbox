function [y,yr,n,s]=tsmav1(tsm,YRS,k1,I1)
% [y,yr,n,s]=function tsmav1(tsm,YRS,k1,I1)
%
% Average of columns in a time-series matrix;  control
% over (1) which variables and years to include, and 
% (2) handling of NaN missing values
%
%*********************** IN ARGS ********************************
%
% tsm (m1 x n1)r time series matrix, m1 years, n1 variables;
%		some elements might be NaN
% YRS (2 x 2)i start, end years of:
%		row 1: tsm
%		row 2: desired output averaged series
% k1 (1 x 1)i option for handling NaN
%		1=NaN's OK for some series; average over however many series
%			for the year are available with valid data
%			disregarding those with NaN
%		2=abort if any NaN's in specified subset of 
%			row and cols marked by YRS and I1 
%		3=set y and s to NaN for any year that not all series specified
%			by I1 have valid data 
% I1 (n1 x 1): mask to use (1) or not use (0) any of the n1 variables
%
%************************ OUT ARGS *********************************
%
% y (my x ny)i averaged series
% yr (my x 1)i year vector for y
% n (my x 1)i  number of series in the average for each year
% s (my x 1)r  standard deviation of chronology vaues for each year;
%		Assigned as NaN if sample size less than 3 chronologies 

a=NaN;

%*** Check inputs

[m1,n1]=size(tsm);
[m10,n10]=size(I1);

if m10~=n1 | n10~=1,
	error('I1 must be cv; length must match col size of tsm')
end
L10=  I1==1 | I1==0;
if sum(L10)~=n1,
	error('I1 must be logical')
end
if  ~any(I1),
	error('At least 1 element of I should be 1')
end


[m10,n10]=size(YRS);
if m10~=2 | n10~=2,
	error('YRS must be 2 x 2')
end
if YRS(2,1)<YRS(1,1) | YRS(2,2)>YRS(1,2),
	error('Period for averaging outside tsm coverage')
end

[m10,n10]=size(k1);
if (m10~=1) | (n10~=1) | (k1 ~=1 & k1 ~=2 & k1~=3),
	error('k1 must be scalar, and be only 1,2 or 3')
end


% Make year pointers
yr1=(YRS(1,1):YRS(1,2))';  % year vector for tsm
L1=yr1>=YRS(2,1) & yr1<=YRS(2,2); % pointer to target years in tsm
yr = yr1(L1);  % year vector for averaged series

% Pull row and column subset to be analyzed
nsers=sum(I1);  % total number of series to be considered for use
X = tsm(L1,I1);

% Compute number of valid series in each year
L2=~isnan(X);   %  pointer to good data
L3=~L2;        %  pointer to NaN's
n =(sum(L2'))'; %  col vector of number of valid data values in each year 
nz=sum(sum(L3)); % total number of NaN's in all years

% Abort if you specified k1==2 and found that row/col subset of data
% contains missing values
if k1==2 & nz~=0; % specified that all series must have valid
	% data in all year, and now see that they do not
	error('NaN data in specified period and k1==2')
end

% For averaging, replace the NaN entries with zeros, compute the
% row sums, then divide by the number of valid variables
% If number of valid values is zero, divide by NaN
X1=X; %operate on a copy of X
if nz>0
	X1(L3)=zeros(nz,1); % X1 has zeros in place of NaNs
end
sum0=sum(n==0); % number of years with zero sample size
LL0=n==0;
sum1=sum(n==1); % number of years with sample size 1
LL1=n==1;
sum2=sum(n>=2);  % number of years with sample size 2 or greater
LL2=n>=2;

%*******************************************
% Compute sample mean. Mean for the year is NaN if sample size zero
nn=n; % work on copy of n that will have NaN in place of zeros
if sum0>0;
	nn(n==0)=a(ones(sum0,1),:);
end
y =  (sum(X1'))' ./ nn; % divide sum of values for each year by
	% either the valid sample size, or NaN if sample size is zero
	% NaN data values do not contribute to sum for the year because
	% those values have been replaced by 0 in X1

%******************************************************
% Compute sample standard deviation of valid values for each year
% Work only on subset of years with 2 or more chronologies
% For years with sample size 1 , set stdev to 0
% For years with sample size 0, set std dev to Nan
% 

% Initialize vector of std devs to NaN
nyrs = length(yr); % number of years of data
s=a(ones(nyrs,1),:); 


% If sample size is 1, set std dev to zero
if sum1>0;
	s(LL1)=zeros(sum1,1);
end


% If sample size is 2 or greater, compute the standard deviation, 
% using the "N-1" denominator in dividing the sums of squares of
% departures
X2=X(n>=2,:); % data for subset of years with sample size >=2

if ~isempty(X2);
	L4=L3(n>=2,:); % recall that L3 is pointer to NaNs
		% L4 is a subset of rows of this pointer
	nz2=sum(sum(L4)); % number of NaN values in all years with
		% sample size >=2
	y2=y(n>=2);  % sample means for the subset years
	n2=n(n>=2); % sample sizes for the subset years

	X2m = y2(:,ones(nsers,1)); % col-dupe the means vector
	XM=X2-X2m; % departures from mean
	XS=XM .* XM;  % squared departures from mean
	XS(L4)=zeros(nz2,1); % replace invalid squared departures with zero
		% so they won't contribute to sum of squares of departures
	xs =   (sum(XS'))'; % cv of sum of squares of departures
	s2 = sqrt(xs ./ (n2-1));  % standard deviation

	s(n>=2)=s2; % store standard devs for this subset of years
		% in vector of std devs for all years
end


% If you set k1==3,  set mean and std dev to NaN for any years
% that have any data missing
L7=n<nsers;
s7=sum(L7);
if k1==3 & s7>0;
	s(L7)=a(ones(s7,1),:);
	y(L7)=a(ones(s7,1),:);
end

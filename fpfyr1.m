function S=fpfyr1(x,yr,fmt1)
%
% Ascii table of an annual time series, one decade per line
%
% Meko 9-8-96
%
%******************* IN ARGS ****************************
%
% x (mx x 1)r   time series
% yr (myr x 1)i  matching year vector for x
% fmt1 (1 x ?)s   format for a data value (say, '%5.2f')
%
%******************* OUT ARG ****************************
%
% S (mS x nS)s   string matrix of output data;note: ms==myr
%

a=NaN;

% Check input
[mx,nx]=size(x);
if mx<2 | nx>1;
	error('x must be col vector')
end
[myr,nyr]=size(yr);
if myr<2 | nyr>1 | myr~=mx
	error('yr must be cv same size as x');
end
if ~isstr(fmt1);
	error('fmt1 must be a string')
end

% Build format for a year and 10 data values
s1=fmt1;
fmt=['%4.0f' s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 '\n'];


% If necessary, tack on leading and trailing NaNs so that the
% time series begins in a "0" year and ends in a "9" year,
% like 1710, 1999
yrgo1=yr(1);
yrsp1=yr(myr);
% How many years to slap in front?
n1=rem(yrgo1,10);
% In back?
n2=10-rem(yrsp1+1,10);
% build the augmented years vectors, if needed
if n1~=0;
	yrb=[(yrgo1-n1):(yrgo1-1)]';
	xb=a(ones(n1,1),:);
	yr=[yrb; yr];
	x= [xb; x];
end
if n2~=0 & n2~=10;
	yre=[(yrsp1+1):(yrsp1+n2)]';
	xe=a(ones(n2,1),:);
	yr=[yr; yre];
	x= [x; xe];	
end	
	

% Shape the matrix to have a year and 10 values per col
mnew=length(x);
ncols=mnew/10; % number of cols in the target matrix X
nrows=mnew/ncols; % number of rows in the target matrix X
X=reshape(x,nrows,ncols);

% Pick every 10th year to put in row 1 of X
w=yr(1):10:yr(length(yr));

% Put actual start and end year in first and last elements of w
%w(1)=yrgo1;
%w(length(w))=yrsp1;

% Slap w as top row on X
X=[w; X];

% Print X to string variable
S=sprintf(fmt,X);

function [X,yr]=sov2tsm(v,YRS,T,J)
% sov2tsm: strung out vector to time series matrix
% CALL: [X,yr]=sov2tsm(v,YRS,T,J)
% Meko 2-20-97
%
%******************** IN **********************************
%
% v (mv x 1)r sov of one or more time series
% YRS(nsers x 3)i start year, end year, and starting row index of
%		each series in v
% T(1 x 2)i <see notes> first last year of desired tsm X 
% J(? x 1)i <see notes>index to rows of YRS specifying which series in 
%		v to include in the tsm.  For example, J=(1:nsers)' would
%		mean include all series, while J=([1 3 4 7])' would pick
%		only those four series.  The series numbers correspond
%		to rows of YRS
%**************** OUT *********************
%
% X(mX x nX)r time series matrix, mX years and nX columns
% yr(mX x 1)i year vector for X
%
%***************** NOTES **************************************
%
% X filled out with NaNs
% X covers period from start of earliest series in v to end of 
%		latest series
% 
% If T is set to [], by default X covers period from earliest to 
%		most recent year in any series in v
% If J is set to [], by default all series in v are included in X
%************************************************************

if nargin~=4;
	error(['Must have 4 input args, has ' int2str(nargin)]);
end

% Size
[mv,nv]=size(v);

% Find earliest an latest year of any series
yrgo = min(YRS(:,1));
yrsp = max(YRS(:,2));
yr=(yrgo:yrsp)';
nyrs=length(yr);

% get number of series
[nsers,nX]=size(YRS);
if nX~=3,
	error('Col size of YRS should be 3');
end



if isempty(T);; % no specified time period for output tsm
	T=[yrgo yrsp];
	kT=0;
else
	kT=1;
end
if isempty(J);
	J=[1:nsers]';
	kJ=0;
else
	kJ=1;
end


% Allocate for X
a=NaN;
X=a(ones(nyrs,1),ones(nsers,1));

% Loop over series
for n=1:nsers;
	yron=YRS(n,1);
	yroff=YRS(n,2);
	% get start and end row index in v
	i1=YRS(n,3);
	i2=i1+(yroff-yron);
	% compute start and end row index in X
	i3 = yron - yrgo +1;
	i4= yroff-yrgo+1;
	% Get the data, and put in X
	x = v(i1:i2);
	X(i3:i4,n)=x;
end

% Trim or expand the rows of the tsm matrix if desired to make it cover
% period specified in T.
if kT==1;  % you specified the time period for the tsm
	% Pull subset of rows of X and T in the specified period
	L1=yr>=T(1,1) & yr<=T(1,2);
	X=X(L1,:);
	yr=yr(L1,:);
	nyrs=length(yr);
	% Leading years (rows) on if needed
	nslap=yr(1)-T(1,1); % need this many years
	if nslap>0;
		yrfront=(T(1,1):(yr(1)-1))';
		Xfront=a(ones(nslap,1),ones(nsers,1));
		X=[Xfront;X];
		yr = [yrfront;yr];
		nyrs=length(yr);
	end
	% Trailing years if needed
	nslap=T(1,2)-yr(nyrs);
	if nslap>0;
		yrback=((yr(nyrs)+1):T(1,2))';
		Xback=a(ones(nslap,1),ones(nsers,1));
		X=[X;Xback];
		yr=[yr; yrback];
		nyrs=length(yr);
	end
end


% Cull selected columns (variables) of X
if kJ==1;
	[mJ,nJ]=size(J);
	if nJ~=1;
		error('J must be cv');
	end
	if any(J>nsers);
		error('An element in J is larger than number of series');
	end
	
	X=X(:,J);
end


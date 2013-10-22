function [Y,yrsy]=seaspt1(X,yrs,msk,ns,nm,sc,k)

% D. Meko 5-7-92

% Given: successive "update" decks of ppt or temp for various stations,
%        all stations covering same years, each record with 13 cols, 
%        year in col 1.

% Gives: Corresponding array of seasonalized data, one col per season


%*******  INPUT ARGS
%
% X (m1 x n1)  monthly climate arrays, stn 1 followed by stn 2, etc
% yrs (1 x 2)  start and end year for the monthly climate series
% 				should be same for each series
% msk (24 x nm) a mask (either 1 or 0) telling which months go into
%    each season.  Each col for a different season.  Row 1 of the
%    24 rows corresponds to january of previous year, row 12 to dec
%    of previous year, row 13 to jan of current year, row 24 to dec of 
%    present year.  For example,
%    [0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0]'  specifies
%    DEC-FEB 
% ns (1 x 1)  number of stations
% nm (1 x 1) number of monthly groupings (seasons)
% sc (1 x 1) scale factor, to be multiplied times resulting seasonalized
%    series to conver units if necessary.  For example, if X are in 
%    hundredths of inches and you want seasonalized in inches, set
%    sc = 0.01.  If no change desired, use sc = 1.0.
% k  (1 x 2)  options
%    k(1) =1 seasonalized is sum of monthly
%         =2 seas         is ave of monthly -- for temperature
%    k(2) =1 if any of the seasonal masks cross the year boundary
%         =0 if not


%************     OUTPUT
%
% Y (m2 x nm)  seasonalized climate array, each col a season, each
%		row a year.  Year keyed to calender year of highest row of
%		msk containing   1.  For example, if yrs says that monthly data
%     covers 1901-1978, and if the mask has any 1s in rows before
%     row 12, the seasonalized array begins in yrs(1)+1.  If all
%     nonzero elements of msk are confined to rows 13-24, then
%     seasonalized array begins in yrs(1).  It is possible to specify
%     a previous fall and a current fall.  Say fall is sept-nov.  
%     If msk ([9 10 11],?) =1, previous fall
%     If msk ([21 22 23],?)=1, current fall
% yrout (1 x 2) begin and end years of seasonalized array Y


%*************  NOTES
%
% You typically should have run gutsrch1.f to get the file X
%
% Good idea to save yrout with Y, and to have some text file with
% description of contents of Y


if nargin~=7; error('NEED 7 INPUT ARGS');  end;
if nargout~=2; error('NEED 2 OUTPUT ARGS');  end;

[m1,n1]=size(X);
nyrs2=yrs(2)-yrs(1)+1;  % number of years of monthly climate input


%***************  CONSISTENCY CHECKS
%
L1 = [nyrs2*ns~=m1  n1~=13];
if any(L1);  error('# ROWS, COLS, OF X INCONSIST WITH ns and yrs'); end;

%  Check that each station has same first and last data year

yrs1=yrs(1);  yrs2=yrs(2);
fst=yrs1(ones(ns,1),:)
lst=yrs2(ones(ns,1),:)

I2=(nyrs2:nyrs2:m1)'
I1=I2 - (nyrs2-1)

L2= fst==X(I1,1);
L3= lst == X(I2,1);
pause

L4=[sum(L2)~=ns  sum(L3)~=ns]
if any(L4); error('NOT ALL STATIONS HAVE SAME START, END YR'); end;

%*******  DIFFERENT BEGINNING YEARS DEPENDING ON k(2)

if k(2) ==1 ;  % at least one msk value on "previous" year
	yrgo=yrs(1)+1;
	I3= (1:nyrs2-1)';
	I4= I3+1;
	nyrs = nyrs2-1;
elseif k(2) == 0;
	yrgo=yrs(1);
	I3=(1:nyrs2); 
	I4 = I3;
	nyrs=nyrs2;
else
	error('k(2) must = 1 or 0');
end

Y = zeros(nyrs*ns,nm);  % preallocate


for i=1:nm;  % loop for each season
	size(msk)
	i
	
	masks =(msk(:,i))';  % month mask for this season
	sm = sum(masks); % number of months in season
	for j=1:ns;  % loop for each station
		I5 = (j-1) * nyrs + (1:nyrs)';
		inc = nyrs2* (j-1);
		X1 = [X(I3+inc,2:13)   X(I4+inc,2:13)];
		X2 = X1(:,masks);
		X3 = (sum(X2'))';
		if k(1) ==2
			X3 = X3 / sm;
		else
		end
		Y(I5,i) = X3 * sc;
	end
end
yrsy = [yrs(1)+k(2)  yrs(2)];


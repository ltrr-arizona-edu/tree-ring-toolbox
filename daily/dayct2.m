function dayct2
% dayct2:  count number of gridpoint or region rainy days each year in specified seasonal window
% dayct2;
% Last revised:  11-26-01 to accomodate optionally mean-adjusted files
%
% Previously ran gridday1.m to get storage file (e.g., day1b.mat) with logical rainy-or-not daily 
% tsm for gridpoints, and a daily tsm of fraction of gridpoints rainy.  Now want to summarize data
% over season (e.g., summer).  Question is how many rainy days each year in that season at each gridpoint.
% Another question is how many regionally rainy days each year, based on a threshold fraction of 
% gridpoints rainy. And, total precip each season.
%
%*** INPUT
%
% No input args
% Prompted for infile produced by gridday1.m
% Prompted for name of outfile to hold results
%
%
%*** OUTPUT
%
% No output args
% Outfile of results has data on number of rainy days each year in the seasonal window, with
%   associated data.  See vlist stored with outfile for definitions
%
%*** REFERENCE -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Input daily tsm should be in the usual form. Year, month, day, "366"-day as cols 1-4.  Time series for
% gridpoints, etc. as remaining cols.
%
% Output time series matrix is an annual time series matrix.  Col 1 is the year, defined as the year of the 
% ending day of the specified period.  So if the season crosses from Dec into Jan, year is that of Jan.
% Example:  [241 31] is previous year's Sept 1 (about) to Jan 31
% Remaining cols are the annual values of number of rainy days, or seasonal total precip
%
%
%***********************  NOTES  ***********************************
%
% Template: dayct1.m


% Previously name rainy.m when developed for Connie Woodhouse rainy day analysis
%
% Daily ppt matrix X should contain 366 days for each year, with possibly
% fewer days in the first and last year.  It is assumed that all missing 
% data have been filled in beforehand.  Day 60 (Feb 29) on non-leap years
% should have a missing-value code.  It makes no difference to this function
% what the code is, because the code is replaced with zeros for the count
% of number of rainy days and the total rainfall over the day window.
%
% Be careful to consider the time coverage by valid daily data in choosing
% the analysis period YRS(2,:).  For example, if valid daily data start in
% Jan 5, 1900 and end in Nov 28, 1993, the combination 
%  endday = [1 366], YRS = [1900 1993; 1900 1993] 
% gives bogus values for years 1900 and 1993.  If you want to have a day
% window 1-366 for this example, you can at best get an annual series
% for 1901-1993 (i.e., YRS=[1900 1993;  1901 1992]).
%
%
%******************  OUT ARGS ***********************************
%
% Y1 (mY1 x nY1)   number of rainy days in day window (enddays)
%		at each station in each year of the specified analysis period 
%		col 1:  year
%		col 2-nY1:  values for each of nY1-1 stations
% 		rows:  each row is a year
% Y2 (mY2 x nY2)   like Y1, except total rainfall, in inches
%
%
%

% PLAN
%
%-- Prompt for specs
%-- Get input data file
%-- Flesh out daily matrix to begin Jan 1, end Dec 31
%-- Compute year information for the output day-windowed series
%-- Compute long-term daily mean of fraction of gridpoints rainy

%-- Compute n-rainy series for regional fraction-rainy data
%-- Compute n-rainy series for individual gridpoints


%-- Prepare output

%---- Prompt for whether using non-adjusted data or means-adjusted
kadj = menu('Choose',...
    'Unadjusted daily ppt as input',...
    'Adjusted (means-adjusted) daily ppt as input');
if kadj==1;
    meanadj='No';
else;
    meanadj='Yes';
end;



%---- Prompt for start and end day of season
prompt={'Enter start day','Enter end day'};
def={'183','244'};
dlgTitle='Starting and ending days of season';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
endday=[str2num(answer{1})  str2num(answer{2})];
clear def dlgTitle lineNo answer 

%---- Prompt for fraction of gridpoints wet to define regionally wet day
prompt={'Enter fraction'};
def={'.50'};
dlgTitle='Threshold fraction of gridpoints wet defining regionally wet day';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
fract=[str2num(answer{1}) ];
clear def dlgTitle lineNo answer 


%-- Get input data file
switch meanadj;
case 'No';
    [file1,path1]=uigetfile('day?b.mat','Infile with gridpoint daily tsm of wet or not');
case 'Yes';
    [file1,path1]=uigetfile('day?bx.mat','Infile with gridpoint daily tsm of wet or not');
otherwise;
end;


pf1=[path1 file1];
eval(['load ' pf1 ' G H pcrit U;']);
% G is daily tsm of wet or not for each gridpoint;  gridpoint logical variables in cols 5-on
% H is cv of daily fraction of gridpoint wet;  same row size as 
% U is cv of daily total precip for the region, computed from griddpoint data by griddata1

%-- prompt for name of output file
switch meanadj;
case 'No';
    deffn = ['day' file1(4) 'c.mat'];
case 'Yes';
    deffn = ['day' file1(4) 'cx.mat'];
otherwise;
end;

[file2,path2]=uiputfile(deffn,'Outfile to hold summary number of rainy day data');
pf2=[path2 file2];


%-- Flesh out daily matrices G and U to begin Jan 1, end Dec 31

[mG,nG]=size(G);
daygo=G(1,4);
yrgo = G(1,1);
daysp=G(mG,4);
yrsp = G(mG,1);
mss = NaN;
if daygo~=1 ; % first row of G is not for Jan 1
	daycv = (1:(daygo-1))';  % cv of days to splice at start of X
	nadd = length(daycv);
	yrcv = repmat(yrgo,ndd,1);  % cv of duped first year
	MSgo = repmat(NaN,nadd,(nG-4));
    Ggo = [yrcv repmat(NaN,nadd,2) daycv MSgo];
	G=[Ggo; G];
    H=[repmat(NaN,nadd,1) ; H];
    U=[repmat(NaN,nadd,1) ;  U];
        
end
if daysp~=366 ; % last row of G is not for Dec 31
	daycv = ((daysp+1):366)';  % cv of days to splice at end of G
	nadd = length(daycv);
	yrcv = repmat(yrsp,nadd,1);  % cv of duped last year
	MSsp = repmat(mxx,nadd,(nG-4));
	Gsp = [yrcv daycv MSsp];
	G=[G; Gsp];
    H=[H;  repmat(NaN,nadd,1)];
    U=[U ; repmat(NaN,nadd,1) ];
    
    
end

% Check data consistency
% Compatible row size of X and YRS?
[mG,nG]=size(G);
if rem(mG,366)~=0;
    error('row size of G should be even multiple of 366');
end;


%-- Compute year information for the output day-windowed series

year = G(:,1); % cv of years for input daily tsm
% Possible first year of day-windowed output, and first needed year
% of input data depend on whether day window crosses the calendar year
% boundary.
k1 = 0;  % indicator for day window crossing year boundary
if  endday(1) < endday(2);  % day window does not cross calendar year boundary
	first = year(1);  % will need data from this year to form first year's grouping
	goposs = year(1);  % first possible year for day-window data
else;  % crosses year boundary
	first = year(1) - 1;  % will need preceding years data
	goposs = year(1) + 1;  % earliest possible year for day-windowed series
	k1 =1 ;   % flag indicating that day window crosses year boundary
end;
spposs=max(year); % last possible year of day-windowed series
nyrs = spposs-goposs+1;  % # of years of output series
yr = (goposs:spposs)';  % cv of output years
nstns = nG -4;  % number of gridpoint series


%-- Compute long-term daily mean of fraction of gridpoints rainy

Hmn=repmat(NaN,366,2);
Hstd=repmat(NaN,366,2);
Hmn(:,1)=(1:366)';
Hstd(:,1)=(1:366)';
for i=1:366
   L5=G(:,4)==i; % Select this day; H is same row-arrangement as G
   HH1=H(L5);
   hmn=nanmean(HH1);
   Hmn(i,2)=hmn;
   hstd=nanstd(HH1);
   Hstd(i,2)=hstd;
end


%-- Compute long-term daily mean of total regional pcp 

Umn=repmat(NaN,366,2);
Ustd=repmat(NaN,366,2);
Umn(:,1)=(1:366)';
Ustd(:,1)=(1:366)';
for i=1:366
   L5=G(:,4)==i; % Select this day; H is same row-arrangement as G
   UU1=U(L5);
   umn=nanmean(UU1);
   Umn(i,2)=umn;
   ustd=nanstd(UU1);
   Ustd(i,2)=ustd;
end


%-- Compute annual time series of fraction of gridpoints rainy in time window


%---- Compute number of days in day window, and form pointer to days
if k1==0;  % not cross year boundary
    L2 = G(:,4) >= endday(1) & G(:,4)<= endday(2);
    L6 = Hmn(:,1) >= endday(1) & Hmn(:,1) <= endday(2);
    ndays = endday(2) - endday(1) +1;  % # days in day window
else
    L2=(G(:,4)>=endday(1) & G(:,4)<=366) | (G(:,4)>=1 & G(:,4)<= endday(2));
    L6=(Hmn(:,1)>=endday(1) & Hmn(:,1)<=366)|(Hmn(:,1)>=1 & Hmn(:,1)<= endday(2));
    ndays = (366 - endday(1)+1) + endday(2);
end

% Make pointer to years of G needed for analysis period;  depends on whether season crosses Dec 31, via first
L1 = G(:,1) >= first & G(:,1) <= spposs;

% allocate for time series & initialize
% F: col 1= year
%   col 2 = number of days with fraction of gridpoints wet >= fract
%   col 3 = average daily fraction of gridpoints wet
F=[yr repmat(NaN,nyrs,3)]; 

%  Get needed rows of H , Hmn and Hstd,etc
H1 = H(L1 & L2,:);
Hmn1=Hmn(L6,2);
Hstd1=Hstd(L6,2);

U1 = U(L1 & L2,:);
Umn1=Umn(L6,2);
Ustd1=Ustd(L6,2);

% Truncate leading and trailing if window crosses year boundary;
if k1==1;
	H1(1:endday(2),:)= [];
	[mH1,nH1]=size(H1);
	nout = 366 - endday(1) + 1;  % # trailing days to drop off
	nn = mH1 - nout;
	H1=H1(1:nn,:);
    
    U1(1:endday(2),:)= [];
	[mU1,nU1]=size(U1);
	nout = 366 - endday(1) + 1;  % # trailing days to drop off
	nn = mU1 - nout;
	U1=U1(1:nn,:);
    
    
end
% H1, U1 should now have only the specified days in the day window
[mH1,nH1]=size(H1);
nr2 = ndays * (spposs-goposs + 1);  % expected number of rows in H1
if nr2 ~= mH1
	error('Row size of H1 incorrect')
end


%---- Compute number of valid (non-NaN) days each year in the endday window

Htemp= reshape(H1,ndays,nyrs);
nvalidH = (sum(~isnan(Htemp)))';
nvalidH=[yr nvalidH];


Utemp= reshape(U1,ndays,nyrs);
nvalidU = (sum(~isnan(Utemp)))';
nvalidU=[yr nvalidU];


%---- COMPUTE ANNUAL TIME SERIES OF NUMBER OF REGIONALLY WET DAYS IN SPECIFIED SEASONAL WINDOW

NH =[yr  (sum(Htemp>fract))'];
Z = [yr  (sum(Utemp))'];


%
% Recall a regionally wet day defined as fraction of gridpoints wet >= fract
%   where a gridpoint/day is defined as wet if the median station P of stations in locus wet, and
%   and a station is defined as wet if daily P>pcrit, as stored in input file day?b.mat
%   

%--- OPERATE ON INDIVIDUAL GRIDPOINTS
%
% Recall that G is a logical daily tsm of rainy days at individual gridpoints. Cols 5-on hold
% the data for each gridpoint.  See pf1 input file variable Gxy(Ireg,:) for the gridpoints, and
% Gstn{} for the stations for each gridpoint

%-- For each gridpoint, compute long-term fract of years wet, for each day of year

GF=repmat(NaN,366,nstns+1);
GN=repmat(NaN,366,nstns+1);
GF(:,1)=(1:366)';
GN(:,1)=(1:366)';
for i=1:366
   L5=G(:,4)==i; % Select this day
   GG1=G(L5,5:(nstns+4));
   GF(i,(2:(nstns+1))) =  nanmean(GG1);
   % Compute number of valid (non-NaN) days, i.e., how many years the mean for that day is based on
   GN (i,(2:nstns+1)) =    sum(~isnan(GG1));
end;



%  The annual tsm of gridpoint number of rainy days in window

%  Get needed rows of G 
G1 = G(L1 & L2,5:(nstns+4));


% Truncate leading and trailing if window crosses year boundary;
if k1==1;
	G1(1:endday(2),:)= [];
	[mG1,nG1]=size(G1);
	nout = 366 - endday(1) + 1;  % # trailing days to drop off
	nn = mG1 - nout;
	G1=G1(1:nn,:);
end
% G1 should now have only the specified days in the day window
[mG1,nG1]=size(G1);
nr2 = ndays * (spposs-goposs + 1);  % expected number of rows in G1
if nr2 ~= mG1
	error('Row size of G1 incorrect')
end

% Allocate & initialize for annual tsm of number of days wet at each gridpoint
% K: col 1= year
%   col 2 = number of days wet at gridpoint 1, in order as Gxy(Ireg,:)
%   col 3 = number ... gridpoint 2 ...
K=[yr repmat(NaN,nyrs,nstns)]; 
KN = [yr repmat(NaN,nyrs,nstns)]; % number of valid days (might be less than ndays) elements of K computed on


for n = 1:nstns ;   % loop over gridpoints
    x2 = G1(:,n);  %  get a gridpoint's daily data indicator of wet or not
    x2 = reshape(x2,ndays,nyrs);  % make into a matrix convenient for summing
    
    K(:,n+1)=(nansum(x2))';
    KN(:,n+1)=  (sum(~isnan(x2)))';
end;



% SAVE OUTPUT

vlist = ['Produced by dayct2.m  on ' pf1 ' with'];
vlist=char(vlist,['   pcrit = ' num2str(pcrit) ' : pre-set threshold (in) for wet day at each station']);
vlist=char(vlist,['   endday = ' int2str(endday) ' : start and end day of season']);
vlist=char(vlist,['   fract = ' num2str(fract) ' : critical fraction of gridpoints required wet']);
vlist=char(vlist,['Hmn long-term mean daily fraction of gridpoints wet']);
vlist=char(vlist,['Hstd ... standard dev of daily ....']);
vlist=char(vlist,['nvalidH= number of non-NaN days in seasonal window each year for input time series H']);
vlist=char(vlist,['NH = number of regionally rainy days in seasonal window each year']);
vlist=char(vlist,['GF= fraction of years in which each day of year was wet at each gridpoint (column)']);
vlist=char(vlist,['GN= sample size (number of valid years) on which fractions in GF computed']);
vlist=char(vlist,['K = annual tsm of number of days wet in seasonal window, by gridpoint']);
vlist=char(vlist,['KN = sample size (number of non-NaN days) on which sums in K computed']);
vlist=char(vlist,['Z = annual tsm of total regional precip in seasonal window']);
vlist=char(vlist,['Date produced = ' date]);


% vlist=char(vlist,
set1 = ' vlist endday pcrit fract Hmn Hstd nvalidH NH GF GN K KN Z';

eval(['save ' pf2  set1 ';'])


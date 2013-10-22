function dayct3
% dayct3:  tsm total regional P and 3 rainy days in a variable time window
% dayct3;
% Last revised:  11-26-01
%
% Previously ran gridday1.m to get storage file (e.g., day1b.mat) with logical rainy-or-not daily 
% tsm for gridpoints, and a daily tsm of fraction of gridpoints rainy.  This aimed at later finding out 
% if strength of signal varies by 1) starting day of window, and 2) number of pentads (5 days) in window
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
% Output is two time series matrices.  These are annual matrices.  Col 1 is the year, defined as the year of the 
% ending day of the specified period.  So if the season crosses from Dec into Jan, year is that of Jan.
% Example:  [241 31] is previous year's Sept 1 (about) to Jan 31.  Col 2 is summary for first grouping, col 2 for 
% second, etc.
% There are npentad groupings
% Bday1, Eday1 are cv's of start and end day for the early-period windows
% Bdat2, Edat2 ... later-period windows
% A is the annual seasonal window tsm for total precip
% B is ... for number of rainy days
%
%***********************  NOTES  ***********************************
%
% Template: dayct2.m
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

%-- Prompt for specs


%---- Prompt for whether using non-adjusted data or means-adjusted
kadj = menu('Choose',...
    'Unadjusted daily ppt as input',...
    'Adjusted (means-adjusted) daily ppt as input');
if kadj==1;
    meanadj='No';
else;
    meanadj='Yes';
end;


%---- Prompt for start days -- one earlier than the other, and maximum number of pentads (5-day periods)
prompt={'Enter first candidate start day','Enter second candidate  start day','Number of pentads'};
def={'173','183','18'};
dlgTitle='Candidate starting days of season, and # pentad';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
dayon=[str2num(answer{1})  str2num(answer{2})];
npentad = str2num(answer{3});
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
    deffn = ['day' file1(4) 'd.mat'];
case 'Yes';
    deffn = ['day' file1(4) 'dx.mat'];
otherwise;
end;
[file2,path2]=uiputfile(deffn,'Outfile to hold summary number of rainy day data');
pf2=[path2 file2];


%--- Compute starting and ending day for each of the npentad pentads
Bday1 = repmat (dayon(1),npentad,1);
Bday2 = repmat (dayon(2),npentad,1);
Eday1 = ((1:18)*5)' + Bday1 - 1;
Eday2 = ((1:18)*5)' + Bday2 - 1;


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
if all(Bday1 < Eday1) & all(Bday2<Eday2);  % no day window crosses calendar year boundary
	first = year(1);  % will need data from this year to form first year's grouping
	goposs = year(1);  % first possible year for day-window data
else;  % crosses year boundary
    error('I did not allow for window crossing year boundary in this code');
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

% Size matrices for output time series
A1 = repmat(NaN,nyrs,npentad+1);
B1 = repmat(NaN,nyrs,npentad+1);
A2 = repmat(NaN,nyrs,npentad+1);
B2 = repmat(NaN,nyrs,npentad+1);
A1(:,1)=yr;
B1(:,1)=yr;
A2(:,1)=yr;
B2(:,1)=yr;


%-- Compute annual time series of fraction of gridpoints rainy in time window


for nn = 1:2;  % loop over the earlier and later periods
    if nn==1; 
        Bday=Bday1;
        Eday=Eday1;;
    else;
        Bday=Bday2;
        Eday=Eday2;
    end;
      
    for n = 1:npentad; % loop over pentads
        
        endday(1)=Bday(n); % start day of window
        endday(2)=Eday(n); % last day of window
                
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
        
        
        %---- COMPUTE ANNUAL TIME SERIES OF NUMBER OF REGIONALLY WET DAYS AND TOTAL REGIONAL RAINFALL
        
        NH =[yr  (sum(Htemp>fract))']; % # rainy days
        Z = [yr  (sum(Utemp))']; % total precip
        
        if nn==1; 
            A1(:,n+1)=Z(:,2);
            B1(:,n+1)=NH(:,2);
        else;
            A2(:,n+1)=Z(:,2);
            B2(:,n+1)=NH(:,2);
        end;
            
        
        %
        % Recall a regionally wet day defined as fraction of gridpoints wet >= fract
        %   where a gridpoint/day is defined as wet if the median station P of stations in locus wet, and
        %   and a station is defined as wet if daily P>pcrit, as stored in input file day?b.mat
        %   
    
        
    end;  % for n = 1:npentad; % loop over pentads
end;  % for nn = 1:2;  % loop over the earlier and later periods


% SAVE OUTPUT
vlist = ['Produced by dayct3.m  on ' pf1 ' with'];
vlist=char(vlist,['   pcrit = ' num2str(pcrit) ' : pre-set threshold (in) for wet day at each station']);
vlist=char(vlist,['   endday = ' int2str(endday) ' : start and end day of season']);
vlist=char(vlist,['   fract = ' num2str(fract) ' : critical fraction of gridpoints required wet']);
vlist=char(vlist,['   dayon  = ' num2str(dayon) ' : start day of the candidate periods']);
vlist=char(vlist,['   npentad = ' num2str(npentad) ':  number of pentads in widest widow']); 
vlist=char(vlist,['   Bday1, Eday1 = start and end day of windows for earlier start day']);
vlist=char(vlist,['   Bday2, Eday2 = start and end day of windows for later start day']);
vlist=char(vlist,['A1 = annual tsm, earlier start day, of total region P for each time window (col)']);
vlist=char(vlist,['A2 = annual tsm, later start day, of total region P for each time window (col)']);
vlist=char(vlist,['B1 = annual # rainy days, earlier start day, in  each time window (col)']);
vlist=char(vlist,['B2 = annual # rainy days, later start day,  in  each time window (col)']);
vlist=char(vlist,['Date produced = ' date]);


% vlist=char(vlist,
set1 = ' vlist endday pcrit fract dayon npentad Bday1 Bday2 Eday1 Eday2 ';
set2 = ' A1 A2 B1 B2';

eval(['save ' pf2  set1  set2 ';'])


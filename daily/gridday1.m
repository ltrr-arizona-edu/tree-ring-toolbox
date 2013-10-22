function gridday1
% gridday1: grid station daily pcp matrix to wet or dry at gridpoints, for later use in computing rainy day frequency
% griday1;
% Last revised 11-27-01
% 
% Have daily pcp matrix for stations in a region.  Want ultimately to get a time series of number of rainy days each year in 
% region.  Will do it via gridpoints. Must first compute gridpoints, then assign each day/gridpoint as wet or not using
% logical variable.  Also (revision) want matrix of daily gridpoint precip
%
%*** INPUT
%
% No input args
% Prompted for 
% 1) file of daily stn precip (e.g., in day1a.mat), at all stations in region.  This was produced by
%   pass 2 of day2fle1.m.
% 2) specs for the grid
%
%
%*** OUTPUT
%
% No output args
% Prompted for:
% 1) file to store logical "wet day" matrix for gridpoints in region
%
%
%*** REFERENCE -- NONE
%*** UW FUNCTIONS CALLED 
%
% gridset
% gcdist 
%
%*** TOOLBOXES NEEDED -- MAPPING
%
%*** NOTES
%
% Plan
%--- Prompt for threshold & gridpoint info.
%--- Get the daily station data 
%--- Conpute the gridpoints, and distances, stations to gridpoints
%-- Associate a nearest gridpoint with each station
%-- Loop over gridpoint
%   Pull submatrix of station data
%   

%-- Save output

%---- Prompt for whether using non-adjusted data or means-adjusted
kadj = menu('Choose',...
    'Unadjusted daily ppt as input',...
    'Adjusted (means-adjusted) daily ppt as input');
if kadj==1;
    meanadj='No';
else;
    meanadj='Yes';
end;


%--- Prompt for threshold & gridpoint info.

prompt={'Enter threshold'};
def={'.05'};
dlgTitle='Threshold (inches of pcp) for a wet day';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
pcrit= str2num(answer{1}); % critical threshold for wetness (in)
pcrit100=pcrit*100;  % mult by 100 because input data in hundredths of inches
clear def dlgTitle lineNo answer 


prompt={'Enter gridpoint spacing, in decimal degrees'};
def={'.25'};
dlgTitle='Desired gridpoint spacing';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
delx = str2num(answer{1});
clear def dlgTitle lineNo answer 


%--- Prompt for types of missing data to accept for rainy days series

kmen1=menu('Choose which type of estimated data to accept for rainy day series',...
    'Real data only (code=0)',...
    'All types of esimated data',...
    'All estimated values, execept by <climatology>');
if kmen1==1; 
    strhow='Rainy day computation using real data only (code==0 daily)';
elseif kmen1==2;
    strhow='Rainy day computation real data and any type of estimate (codes 1-8)';
elseif kmen1==3;
    strhow ='Rainy day computation uses real data plus all estimates except climatology (code==7)';
end;



%--- Prompt for types of missing data to accept for total P amount series

kmen2=menu('Choose which type of estimated data to accept for total precip series',...
    'Real data only (code=0)',...
    'All types of esimated data',...
    'All estimated values, execept by <climatology>');
if kmen2==1; 
    strtmp='Total pcp computation using real data only (code==0 daily)';
elseif kmen2==2;
    strtmp='Total pcp computation real data and any type of estimate (codes 1-8)';
elseif kmen2==3;
    strtmp ='Total pcp computation uses real data plus all estimates except climatology (code==7)';
end;
strhow = char(strhow,strtmp);



%--- Prompt for mean or median in total-precip computations

kmen3=menu('Choose how gridpoint total pcp computed from stations, and region from gridpoints',...
    'nanmean',...
    'nanmedian');
if kmen3==1; 
    strtmp='Total pcp: nanmean used to compute gridpoint from station, region from gridpoint';
elseif kmen3==2;
    strtmp='Total pcp: nanmedian used to compute gridpoint from station, region from gridpoint';
end;
strhow = char(strhow,strtmp);



%--- Get the daily station data 
switch meanadj;
case 'No';
    [file1,path1]=uigetfile('day?a.mat','Infile with daily station data');
case 'Yes';
    [file1,path1]=uigetfile('day?ax.mat','Infile with daily station data');
otherwise;
end;
pf1=[path1 file1];
eval(['load ' pf1 ' stninf T W;']);
[mT,nT]=size(T);
Porig=T; % use this version for total precip
disp([pf1 ' loaded']);



%--- Prompt for output file
switch meanadj;
case 'No';
   [file2,path2]=uiputfile('day?b.mat','Outfile with gridpoint daily data'); 
case 'Yes';
   [file2,path2]=uiputfile('day?bx.mat','Outfile with gridpoint daily data');
otherwise;
end;
pf2=[path2 file2];




%--- Cull acceptable estimated data

if kmen1==1; % accept only observed data for rainy days
    Tslab = T(:,5:nT);
    Lmask = W(:,5:nT) ~=0;
    Tslab(Lmask)=NaN;
    T(:,5:nT)=Tslab;
    clear Tslab;
    clear Lmask
elseif kmen1==2; % all estimates OK 
    % No action needed; T unchanged
elseif kmen1==3; % all types except climatology estimates OKM
    Tslab = T(:,5:nT);
    Lmask = W(:,5:nT) ==7 | W(:,5:nT) ==8;
    Tslab(Lmask)=NaN;
    T(:,5:nT)=Tslab;
    clear Tslab;
    clear Lmask
end;
    

%--- Cull acceptable estimated data for total P

if kmen2==1; % accept only observed data 
    Pslab = Porig(:,5:nT);
    Lmask = W(:,5:nT) ~=0;
    Pslab(Lmask)=NaN;
    Porig(:,5:nT)=Pslab;
    clear Pslab;
    clear Lmask
elseif kmen2==2; % all estimates OK 
    % No action needed; Porig unchanged
elseif kmen2==3; % all types except climatology estimates OKM
    Pslab = Porig(:,5:nT);
    Lmask = W(:,5:nT) ==7;
    Pslab(Lmask)=NaN;
    Porig(:,5:nT)=Pslab;
    clear Pslab;
    clear Lmask
end;


%--- Conpute the gridpoints and distance from stations to gridpoints

lat1=str2num(stninf(:,35:40));
lon1=str2num(stninf(:,41:48));
Gxy=gridset(lat1,lon1,delx);
ngrid=size(Gxy,1); % number of gridpoints
A=fliplr(Gxy); % puts lon in left col
B=[lon1 lat1];
D =gcdist(B,A); % distance, station to gridpoints; each row for a station



%-- Associate a nearest gridpoint with each station

[s1,i1]=min(D'); % rv, distance to nearest gridpoint (km), and which gridpoint that is (index into X);
% rv i1 is index to rows of Gxy telling nearest gridpoint this station


%-- Find the  active (represented) gridpoints -- those with stations in their search radius (in any day)
i2 = unique(i1);
nlive=length(i2);
if nlive<2;
    error('Fewer than two gridpoints in region with any attached stations');
end;



%-- Loop over gridpoints

% Allocate G; first find number of "live" gridpoints
G=repmat(NaN,mT,length(i2));  %to hold # rainy days at gridpoint
P=repmat(NaN,mT,length(i2));  %to hold toal pcp "at" gridpoint
nstn=zeros(nlive,1);  % number of stations "at" each gridpoint
Gstn=cell(nlive,1); % station info at gridpoint

ipoint=0;
for n = 1:nlive; % loop over live gridpoints
    disp(['On gridpoint ' int2str(n)]);
    yx = Gxy(i2(n),:); % lat long coordinates of stations for  gridpoint
    icol=find(i1==i2(n)); % col index in T, Porig, of stations for gridpoint
    nstn(n)=length((icol)); % number of stations for this gridpoint
    Gstn{n}=stninf(icol,:); % station infor for this gridpoint
    T1=T(:,icol+4); % pull submtx of station data for rainy days
    P1 = Porig(:,icol+4); % submatrix for total precip
    LNaN = isnan(T1);
    P1NaN = isnan(P1);
    T1=T1>pcrit100; % Logical 1 if wet
    if any(any(LNaN)); % Want NaN daily stn data to show as NaN, not 0
        T1(LNaN)=NaN;
    end;
    
    % Compute gridpoint series
    if nstn(n)>1;
        G(:,n)=(nanmedian(T1'))';  % cv, 1 if median of stations in locus wet, 0 otherwise
        if kmen3==1; % mean
            P(:,n)=(nanmean(P1'))';  % cv, mean of daily P of stations in locus 
        elseif kmen3==2; % median
            P(:,n)=(nanmedian(P1'))';  % cv, mean of daily P of stations in locus 
        end;
    else;
        G(:,n)=T1;  % special case of just one station attached to gridpoint
        P(:,n)=P1; 
    end;
    
   
end;



%--  Compute fraction of valid (non-NaN) gridpoints wet each day

N =    (sum(~isnan(G')))';  % cv, number of valid gridpoints each day, for rainy days
Ntot = (sum(~isnan(P')))';  % ... for total pcp

% For rainy days, warn if any day has zero valid gridpoints, and is not Dec 29 of a non-leap year
Lcheck=N==0;
if any(Lcheck);
    yr = T(:,1);
    month=T(:,2);
    day=T(:,3);
    L1 =   month==2  & day==29  &   ~leapyr(yr);
    L2 = Lcheck & ~L1;
    if any(L2);
        hwarn1=warndlg('Some day of G-- and not Feb 29 of a nonleap year -- has no valid gridpoints','Warn1');
    end;
    N(Lcheck)=NaN;
end;

H =   (nansum(G'))'  ./  N;  % fraction of live gridpoints wet each day


% For total precip, warn if any day has zero valid gridpoint s and ius not Dec 29 of a non-leap year
Lcheck=Ntot==0;
if any(Lcheck);
    yr = T(:,1);
    month=T(:,2);
    day=T(:,3);
    L1 =   month==2  & day==29  &   ~leapyr(yr);
    L2 = Lcheck & ~L1;
    if any(L2);
        hwarn1=warndlg('Some day of P-- and not Feb 29 of a nonleap year -- has no valid gridpoints','Warn1');
    end;
    Ntot(Lcheck)=NaN;
end;

% Compute regional total precip from gridpoint
if kmen3==1;
    U =   (nanmean(P'))' ;  % mean precip (averaged over active gridpoints with data that day
elseif kmen3==2;
    U =   (nanmedian(P'))' ;  % median precip (of active gridpoints with data that day
end;



T=T(:,1:4);
pack;
G=[T(:,1:4) G ];

Igrid=i2; % active gridpoints (rows of Gxy)


%---- Pull submatrix of station data
  


vlist = ['Produced by gridday1.m  on ' pf1 ' with'];
vlist=char(vlist,['   delx = ' num2str(delx)]);
vlist=char(vlist,['   pcrit = ' num2str(pcrit) ' in']);
vlist=char(vlist,strhow);
vlist=char(vlist,'Gxy lats, lons of gridpoints in region');
vlist=char(vlist,'Gstn -- text info, stations at each live gridpoint');
vlist=char(vlist,'D distance (km) station to gridpoints; 1 row per stn');
vlist=char(vlist,'nstn  number of stations for each gridpoint');
vlist=char(vlist,'Igrid index to rows of Gxy identifying live gridpoints -- data  in cols 5-on of G');
vlist=char(vlist,'G daily tsm of logical wet-or-not for each live gridpoint ; cols 5-on for gridpoints');
vlist=char(vlist,'delx gridpoint spacing (dec deg lat or lon)');
vlist=char(vlist,'pcrit (in) : threshold for a wet day a station'); 
vlist=char(vlist,'N  number of valid (non-NaN) gridpoints each day for rainy day computation');
vlist=char(vlist,'Ntot  number of valid (non-NaN) gridpoints each day for total P computation');
vlist=char(vlist,'H  daily tsm of fraction of live gridpoints wet');
vlist=char(vlist,'P  daily tsm of amount of precip AT gridpoints; mean of station data in locus; same size as G');
vlist=char(vlist,'U  daily tsm of average of gridpoint precip (P) over gridpoints with data ');
vlist=char(vlist,['Date produced = ' date]);

% vlist=char(vlist,


set1=' vlist strhow Gxy G Igrid D nstn Gstn delx pcrit N Ntot H P U';
eval(['save ' pf2  set1 ';']);

disp('here');









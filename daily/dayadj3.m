function dayadj3
% dayadj3: adjust Eischeid daily station precip to reduce bias in regional series due to diffs in means
% dayadj3;
% Last revised 11-1-01
%
% Adjust Eischeid daily station precip to reduce bias in regional series due to differences in means. Written
% for NOAA daily precip analysis to improve stability of regional daily precipitation estimates.  Concerned
% that time series could be biased by variable makeup of stations in region over time.  

%*** NOTES
%
% This version uses Gaussian smoothing of daily data, and a secondary adjustment for reference-period means
% of individual stations for anomalous dryness or wetness of the coverage period
%
%--- STEPS
%
%
% Load daily storage matrix
%
% Initially settings: reference period;  minimum allowable number of years of data in reference period for 
% master stations
%
% Reference period analysis
%   Select master stations
%   Compute regional daily means with master-station dat
%   Smooth (Gaussian + spline) regional daily means
%   Compute and smooth indiv station daily means and compute ratios of region to station
%
% Full-period analysis
%   Scale all daily data
%
% Qualtity control


% Hard Code
nday=[31 29 31 30 31 30 31 31 30 31 30 31]; % n of days in months


% Verbose vs concise mode
kblab = menu('Choose one',...
    'Just the facts, mam',...
    'Verbose (informative-plot) mode');




%--- LOAD DAILY STORAGE MATRIX

[file1,path1]=uigetfile('day?a.mat','Infile of stored daily pcp');
pf1=[path1 file1];
% Check that have duped the day?a.mat as ax 
file2  = [strtok(file1,'.') 'x.mat'];
pf2 =[path1 file2];
if ~(exist(pf2,'file')==2);
    error(['Hey, first must dupe-copy ' pf1 ' as ' pf2 ';']);
end;
    
    
eval(['load ' pf1 ' T W stninf;']);
[mT,nT]=size(T);
nstns=nT-4; % number of stations in region

% Pull off and store 4 leftmost cols 
T1=T(:,1:4);
T(:,1:4)=[];
W1=W(:,1:4);
W(:,1:4)=[];

% Start and end years of input matrix
yrgo=min(T1(:,1));
yrsp=max(T1(:,1));
yr =(yrgo:yrsp)';
nyr = length(yr);

% Convert all estiimated data to NaN
Lbad = W~=0;
T(Lbad)=NaN;


% SET REFERENCE PERIOD AND MINIMIM ALLOWABLE NUMBER OF YEARS OF DATA (ANY DAY) FOR STATION TO BE CONSIDERED MASTER

prompt={'Enter start and end years of reference period',...
        'Enter minimum allow number of years of data in reference period for master station',...
        'Enter period (days) of 0.5 response of guassuian filter',...
        'Enter spline parameter',};
switch(file1(4)); % different settings by region
case '1';
    def={'[1906 1999]','30','19','1E-5'};
    c1=1; % critical precip
case '2';
    def={'[1906 1999]','30','19','1E-5'};
    c1=1;
case '3';
    def={'[1906 1999]','25','19','1E-5'};
    c1=1;
case '4';
    def={'[1920 1999]','25','19','1E-5'};
    c1=1;
case '6';
    def={'[1906 1999]','30','19','1E-5'};
    c1=1;
end;

dlgTitle='Settings';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
yrsref=answer{1};
minref=answer{2};
pgauss = answer{3};
pspln = answer{4};
strack=['Specified start and end years of reference period = ' yrsref];
strack=char(strack,['Specified minimum allow number of years of data in reference period for master station= ' minref]);
strack=char(strack,['Gaussian filter 0.5 response at ' num2str(pgauss) ' days']);
strack=char(strack,['Spline parameter, P=  ' num2str(pspln)]);

yrsref = str2num(answer{1});
minref=str2num(answer{2});
pgauss=str2num(pgauss);
pspln = str2num(pspln);


% Gaussian weights
b=wtsgaus(pgauss);
nb=length(b);
if mod((nb+1),2)~=0;
    error('Gaussian filter has even # weights');
end;

% Store reference-period data
yrgoref = yrsref(1);
yrspref = yrsref(2);
Lref=T1(:,1)>=yrgoref & T1(:,1)<=yrspref;
TB=T(Lref,:);
WB=W(Lref,:);
TB1=T1(Lref,:);
WB1=W1(Lref,:);
yrref =(yrgoref:yrspref)';
nyrref = length(yrref);



%*** REFERENCE PERIOD ANALYSIS


% SELECT MASTER STATIONS -- must have at least minref years of data for each day in ref period

Lmast =logical(zeros(1,nstns)); % initially none are masters
for n = 1:nstns; % loop over stations
    x = TB(:,n); % subset of ref period data for station
    if mod(length(x),366)~=0;
        error(['n of rows of x not evenly divisible by 366 for stn  ' int2str(n)]);
    end;
    ncol = nyrref/366;
    X=reshape(x,366,nyrref);
    X(60,:)=[]; % Feb 29
    nsum1 =  sum(~isnan(X'));
    nsum2= nsum1>minref;
    if all(nsum2);
        Lmast(n)=1;
    end;
end;
nmaster = sum(Lmast);
if nmaster<(0.5*nstns);
    strtemp ={'Fewer than half the total stations are master stations',...
            ['nmaster = ' int2str(nmaster)]};
    uiwait(msgbox(strtemp,'Message','modal'));
    kq = questdlg('Abort the Mission?');
    if strcmp (kq,'Yes') | strcmp(kq,'Cancel');
        error('Mission aborted');
    end;
end;
clear strtemp kq;;
 


% COMPUTE REF PERIOD MEANS 

M1 = repmat(NaN,366,nstns); % means
N1 = repmat(NaN,366,nstns); % sample size
for n = 1:nstns; % loop over stations
    x = TB(:,n); % subset of ref period data for station
    ncol = nyrref/366;
    X=reshape(x,366,nyrref);
    X=X';
    L1=~isnan(X);
    M1(:,n) = (nanmean(X))';
    N1(:,n)= (sum(L1))';
    % Stablilize Feb 29
    x60=[X(:,59); X(:,60); X(:,61)];
    M1(60,n)=nanmean(x60);
end;


% MASTER REGIONAL MEAN -- ys , with day vector ts

Tmast = TB(:,Lmast); % use only the master stations 
ncol = size(Tmast,1)/366;
Tmast = reshape(Tmast,366,(ncol*nmaster));
M2 = (nanmean(Tmast'))'; % regional master daily mean precip, reference period
t=(1:366)';
ttt=(1:(3*366))';
% Gaussian filter 
MMM=[M2 ; M2 ; M2];
[y,ty]=filter1(MMM,ttt,b,1);
y(y<c1)=NaN; % subst NaN for any regional guassian smoothed ppt below threshold
Lf=~isnan(y);
ymaster = y(367:(2*366));
tmaster=(1:366)';
Lmaster = ~isnan(ymaster);
% Spline smooth gaussian smoothed master
% ys = csaps(ty(Lf),y(Lf),pspln,ty);
% ys = ys(367:(2*366));
% ts = ty(367:(2*366));
% plot(t,M2,ty,y,ts,ys);



% COMPUTE AND SMOOTH INDIVIDUAL STATION DAILY MEANS AND COMPUTE RATIOS OF REGION TO STATION (366 X nstns)

R = repmat(NaN,366,nstns); % daily adjustment ratios

for n =1:nstns;
    Tkey = TB(:,n); % store key station data
    ncol = size(Tkey,1)/366;
    Tkey1 = reshape(Tkey,366,ncol);
    m3 = (nanmean(Tkey1'))'; % daily stn means
    t=(1:366)';
    ttt=(1:(3*366))';
    % Gaussian filter 
    MMM=[m3; m3; m3];
    [y,ty]=filter1(MMM,ttt,b,1);
    y(y<c1)=NaN;
    Lfine=~isnan(y);
    ystn = y((367:(2*366)));
    tstn = (1:366);
    Lstn = ~isnan(ystn);
    
    % Ratio of guassian smoothed series
    r2 = ymaster ./ ystn;
    Lboth = ~isnan(r2);
    
    % Smooth ratio 
    R2=[r2; r2; r2];
    plot(ttt,R2);
    LR2 = ~isnan(R2);
    R3= csaps(ttt(LR2),R2(LR2),pspln,ttt);
    r3 = R3(367:(2*366));
    %plot(tmaster,r2,tmaster,r3);
    
    
    
    % Spline smooth
%    ysi = csaps(ty(Lfine),y(Lfine),pspln,ty);
%    ysi = ysi(367:(2*366));
%     yabs =    min(ysi(ysi>0)); % smallest positive smoothed values
%     Lsmall=ysi<=0;
%     if any(Lsmall);
%         ysi(Lsmall)=yabs;
%     end;
    
%    tsi = ty(367:(2*366));
    
    if kblab==2; % & n ==50;
        subplot(2,1,1);
        plot(tmaster,ymaster,tstn,ystn);
        grid;
        legend('Region','Station');
        title(stninf(n,:));
        
    end;
    
%     % Compute and plot ratio
%     r = ys ./ ysi;
%     Ldecent =  ys>c1 & ysi>c1;
%     r(~Ldecent)=NaN;
%     rrr=[r; r; r];
%     Ldecent = ~isnan(rrr);
%     % Spline smooth
%     rs = csaps(ttt(Ldecent),rrr(Ldecent),pspln/10,ttt); % note stronger smoothing for the ratios
%     rs = rs(367:(2*366));
%     rabs =    min(rs(rs>0)); % smallest positive smoothed value
%     Lsmall=rs<=0;
%     if any(Lsmall);
%         rs(Lsmall)=rabs;
%     end;
    if any(r3<=0);
        error(['Zero or lower smoothed daily ratio for ' stninf(n,:)]);
    end;
    R(:,n)=r3; % store adjusment ratios
    
    if kblab==2; % & n==50;   
        subplot(2,1,2);
        plot(t,r2,t,r3);
        ylabel ('Ratio');
        xlabel('Day');
        title(stninf(n,:));
        legend('Before spline','Spline-smoothed');
        grid;
        zoom xon;
        pause;
    end;
    
end;



% ADJUST WHOLE TSM

A = repmat(R,nyr,1); % dupe adjustment matrix
T= T .* A;
T=[T1 T];


% SAVE ADJUSTED DATA  

file2  = [strtok(file1,'.') 'x.mat'];
pf2 =[path1 file2];
eval(['save ' pf2 ' T strack -append;']);






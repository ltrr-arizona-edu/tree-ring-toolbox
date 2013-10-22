function PE=pethorn(T,yrc,lat,dayz,hemis)
% pethorn:  potential evapotranspiration by Thornthwaite method
% CALL: PE=pethorn(T,yrc,lat,dayz,hemis);
%
% Meko 5-23-97
%
%****************  IN **************************
%
% T (? x 13)r    matrix of year (col 1) and 12 monthly Temperature values (oF)
% yrc (1 x 2)i   start and end years of 'calibration' period used to 
%    compute heat index
% lat (1 x 1)i   decimal latitude of the station or center of region
% dayz (51 x 12)r Table 6, p. 228, Thornthwaite & Mather -- table of
%    mean possible monthly duration of sunlight in northern hemisphere
% hemis (1 x 1)s hemisphere (N or S)
%    ( so far, have only keyed in thorthwaites daylenth table for nh)
%********************** OUT **********************
%
% PE (? x 13)r matrix of PE (in);  same size as T
%
%**************** NOTES **************************
%
% Global tables needed: daylen -- these must also be global in calling pgm
%
% Sources:  
%  Sellers, 1960, Physical Climatology, 
%        Heat index equations for I and a from p. 171
%  Pelton, King, and Tanner, 1960, An evaluation of the Thornthwaite and Mean
%        temperature methods for determining potential evapotranspiration. 
%        Argonomy Journal, 387-395 --  eqns 3, 4 for computing undadjusted
%        PE in cm/month and for adjusting for deviation from 12-hr day and
%        30-day month
%  Thornthwaite and Mather 1957, Instructions ...
%        Table 5, p. 226 for unadjusted PE when mean temp above 26.5C
%        Table 6, 7, p. 228,229   mean poss monthly duration of sunlight
%       
%

daysmon=[31 28 31 30 31 30 31 31 30 31 30 31]; % number of days in month

anan=NaN;

switch hemis
case 'N';
   % No action needed
case 'S';
   error('Need to key in dayz and store in daylensh before proceeding');
otherwise;
   error('Invalid hemis');
end


%*************************************************************************
% Build Table for unadjusted PE for t greater than 26.5 oC, or 80 oF. Table
% values are in mm/day, and you specify T in deg C in Thornthwaite.  

Thot = [...
      NaN NaN NaN NaN NaN 4.5 4.5 4.6 4.6 4.6 ...
      4.6 4.7 4.7 4.7 4.8 4.8 4.8 4.8 4.9 4.9 ...
      4.9 5.0 5.0 5.0 5.0 5.1 5.1 5.1 5.1 5.2 ...
      5.2 5.2 5.2 5.2 5.3 5.3 5.3 5.3 5.4 5.4 ...
      5.4 5.4 5.4 5.5 5.5 5.5 5.5 5.5 5.6 5.6 ...
      5.6 5.6 5.6 5.6 5.7 5.7 5.7 5.7 5.7 5.8 ...
      5.8 5.8 5.8 5.8 5.8 5.8 5.9 5.9 5.9 5.9 ...
      5.9 5.9 5.9 5.9 6.0 6.0 6.0 6.0 6.0 6.0 ...
      6.0 6.0 6.0 6.0 6.1 6.1 6.1 6.1 6.1 6.1 ...
      6.1 6.1 6.1 6.1 6.1 6.1 6.1 6.1 6.1 6.1 ...
      6.1 6.1 6.2 6.2 6.2 6.2 6.2 6.2 6.2 6.2 ...
      6.2 6.2 6.2 6.2 6.2 6.2 6.2 6.2 6.2 6.2 ...
      6.2 ];
xThot = 26.0:0.1:38.0;


%********************************************************************

% Check T
[mT,ntemp]=size(T);
if ntemp~=13;
   error('T must be 13=col mtx');
end

% Check years for T data
yr1=T(:,1);
Lleap = rem(yr1,4)==0; % pointer to leap years
yrgo=yr1(1);
nyrs1=length(yr1);
yrsp=yr1(nyrs1);
if yrc(1)<yrgo | yrc(2)>yrsp;
   error('T data does not cover specified normals period');
end

% Compute normals of monthly T for period yrc
L1=yr1>=yrc(1) & yr1<=yrc(2);
TN = T(L1,2:13);
Tmean=nanmean(TN);

% Make matrix of mean T in case need to fill in missing monthly data
TMEAN = repmat(Tmean,mT,1); % this in oF

% Convert rv of normals to centigrade for use in sellers  form of heat index eqn
Tmean=(5/9)*(Tmean-32); % in centigrade


%**************************************

% Apply eqn for heat index I in Sellers, P. 171
% Note that heat index is computed on months whos mean T is above freezing,
% and that for purposes of the computation, any monthly normal above 26.5C is
% set to 26.5C
Tmean1=Tmean(Tmean>0);
Tmean1(Tmean1>26.5)=26.5;
I =    sum((Tmean1/5).^ 1.514);

% Apply eqn for the exponent a
a=1E-6* (0.675*I^3 - 77.1*I^2 + 17920*I + 492390); 


%****************** Compute unadjusted PE in cm/month, approp for 30-day month, 12 hr day

% Pull complete temperature matrix and replace any missing values with 
% monthly mormals
F=T(:,2:13); % data in farenheit
Lmiss = isnan(F);
if any(any(Lmiss));
   F(Lmiss) = TMEAN(Lmiss);
end

% Convert data to centigrade
C = (5/9) * (F-32);


% Pointers to special case data
Lwarm = C>=26.5;
Lhot = C>38.0;
Lcold = C<=0;

% Replace any temperatures above 38 degrees  C with 38 degrees.  This because 
% Table 5, p. 226 of T&M dealing with PE under very high T goes up to 
% 38 deg (100.4F) only.  Note -- I cannot imagine a monthly mean temperare 
% above 38C -- maybe a monthly mean maximum.
nhot = sum(sum(Lhot));
if nhot>0;
   C(Lhot)=38;
end


% Compute unadjusted PE for regular data -- units mm/month
PE = 16 * ((10.0 * C / I) .^ a);

% Replace cold-value PE with zero; PE is set to zero for any T 0C or lower
ncold = sum(sum(Lcold));
if ncold>0;
   PE(Lcold) = zeros(ncold,1);
end


% Handle the high temperature cases, those that need table Thot
nwarm=sum(sum(Lwarm));
if nwarm>0;
   Ctoast = C(Lwarm);
   S = zeros(mT,12);
   % get PE from table
   S(Lwarm) = interp1(xThot,Thot,Ctoast); % Values in mm/day
   % Go from mm/day to mm/month
   DM = repmat(daysmon,mT,1);
   S = S .* DM;
   S(Lleap,2) = S(Lleap,2) * 29/28; % one more day in leap year
   PE(Lwarm) = S(Lwarm); % PE in mm/month
end




%*******************************  ADJUST PE FOR SUNSHINE DURATION 


% The table dayz applies to latitudes 0-50 N, for 12 months of year

% Poleward of 50 deg, use table value for 50 deg
% Interpolate the sunshine adjustment factor from a table
dayfact = interp1((0:50)',dayz,min(lat,50)) / 30; % a rv, length 12
% factor to multiply unadjusted PE by.  

% Expand dayfact to a matrix, same row size as PE
D = repmat(dayfact,mT,1);

% Values in col 2 of D (for Feb) should be based on 29 days instead of 28
% days in leap years.  Change those vales by multiplying by ratio 28/29
D(Lleap,2) = D(Lleap,2) * 28/29;

% Adjust the PE
PE = PE .* D; % PE in mm/month


%*************************  CONVERT from mm/month to inches/month, and slap on yr colmn

PE = PE / 25.4;
PE=[yr1 PE];


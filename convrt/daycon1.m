function [YD,Y,L]=daycon1(stnfn,dwant)
% daycon1: convert US summary of day climate data to Julian day matrix and monthly matrix
% CALL: [YD,Y,L]=daycon1(stnfn,dwant);
%
% Meko 5-12-98
%
%************** IN 
%
% stnfn (1 x ?)s   path\filename of file downloaded from ncdc
% dwant (1 x 1)s  1-letter code for US summary of day telling data type wanted
%    A-daily max tmp (nearest deg F)
%    B-daily min tmp (nearest deg F)
%    C-daily pcp (in, to nearest hundredth)
%    D-daily snowfall (in, to nearest hundredth)
%    E-snow depth at obs time
%    F-max rh
%    G-min rh
%    H-peak gust speed
%    i-daily total sunshine
%    j-daily % of possible sunshine
%    k-departure from normal tmp
%    l-ave daily dewpoint tmp
%    M-heating degree days
%    N-cooling degree days
%    O-ave daily station pressure
%    P-ave daily wind speed
% 
%************** OUT 
%
% YD (? x 3)r   daily Julian day data. 
%  col 1-- year
%  col 2-- julian day
%  col 3 -- daily pcp
%
%******************** NOTES
%
% First written to get climate data formatted for use in Great Sand Dunes
% study.  Alamosa had daily data in US summary of day network. I needed
% monthly totals, and also wanted daily data in convenient format for 
% analysis.
%
% For C,D,M,N data types, monthly value is computed as sum over daily values
% For all other variables, monthly value is mean of daily values
%
% Summary of day data is in English units (e.g., inches, deg F)


% Will want to compare string heading monthly values in file with names of months
strmon={'January','February','March','April','May','June','July',...
      'August','September','October','November','December'};


%------------- OPEN FILE FOR READING
fid1 = fopen(stnfn,'r');


%-------------- READ FIRST 2 LINES AND ALLOCATE STORAGE
c = fgetl(fid1); % should be blank line
c = fgetl(fid1); % line with name of station and period of data
yrgo1 = str2num(c(33:36));
mongo1 = str2num(c(37:38));
yrsp1 = str2num(c(40:43));
monsp1 = str2num(c(44:45));
monthend=strmon{monsp1};
nyr = yrsp1-yrgo1+1; % number of years of data

% Daily julian day matrix
m1 = nyr *366; % number of rows in daily julian day matrix
YD = repmat(NaN,m1,3);
a1=yrgo1:yrsp1;
a1 = repmat(a1,366,1);
YD(:,1)= a1(:);% year column
a1 = (1:366)';
YD(:,2) = repmat(a1,nyr,1);

% Monthly matrix
Y = repmat(NaN,nyr,13);
yr = (yrgo1:yrsp1)';
Y(:,1)=yr;

L=[];

%-------------- READ SECOND BLANK LINE
c = fgetl(fid1);

%----------------- FIND OUT DATA TYPE
c = fgetl(fid1);
dtype = c(2);
if ~strcmp(dwant,dtype);
   disp(c);
   error(['Wrong data type: wanted ' dwant ';  got ' dtype]);
end

%-------------- READ BLANK LINE BEFORE MONTH/YEAR LINE
c = fgetl(fid1);


%************  FILL MONTHLY MATRIX *******************


kfirst = 1;  % flag for first month/year
kwh1=1; % while control
while kwh1==1;
   
   % Read month/year
   c=fgetl(fid1);
   monthis = strtok(c);
   c=fliplr(strtok(fliplr(c)));
   yrthis = str2num(c);
   iadd = (yrthis-yrgo1)*366; % increment later needed for row indexing
   
   %%if strcmp(monthis,'September') & yrthis==1993;
    %%  disp(c);
   %%end
   
   
   % Skip next 3 lines
   c=fgetl(fid1); c=fgetl(fid1); c=fgetl(fid1); 
   % should now be positioned at first daily data for this month/year, or
   % at a blank line if no data for this month/year
   
   kwh2 = 1; % while loop control for data reading
   while kwh2==1;
      c=fgetl(fid1);
      if isempty(c); % reached end of daily data for this month/year
         if yrthis==yrsp1 & strcmp(monthis,monthend);
            kwh1=0;
            kwh2=0;
         else;
            kwh2=0;
            c=fgetl(fid1); % read blank line to position at next month/year line
         end
         
      else; % data hit
         
         % Get day and data value
         c=str2num(c);
         day = c(1);
         x=c(2);
         
         % Compute row slot in storage matrix 
         imonth = find(strcmp(monthis,strmon)); % month index (e.g., 1==jan)
         jday=mdy2jdy(imonth,day); % julian day (1 to 366)
         irow = iadd + jday; % row in target matrix
                          
         
         % Put daily data value in slot
         YD(irow,3)=x;
         YD(irow,:);
         
      end
   end
end; % while kwh1==1



%***************  MONTHLY TOTALS

a1 = [repmat(1,1,31) repmat(2,1,29) repmat(3,1,31) repmat(4,1,30)];
a2 = [repmat(5,1,31) repmat(6,1,30) repmat(7,1,31) repmat(8,1,31)];
a3 = [repmat(9,1,30) repmat(10,1,31) repmat(11,1,30) repmat(12,1,31)];
a4=[a1 a2 a3]';
a4 = repmat(a4,nyr,1);
ndays = [31 29 31 30 31 30 31 31 30 31 30 31];

% Leap year flag
L29 = mod(yr,4) == 0;
L29=L29';
L28 = ~L29;


for n = 1:12; % Loop over months
   L1 = a4==n;
   nd = ndays(n); % number of days in month
   y = YD(L1,3);
   Y1 = reshape(y,nd,nyr);
   if strcmp(dtype,'C') | strcmp(dtype,'D') | strcmp(dtype,'M') | strcmp(dtype,'N');
      %for pcp, snowfall, # heating ddays, # cooling degdays, compute sum over days
      if n==2; % february special: sum over 29 days inleap yr, 28 in other years
         ys29 = (sum(Y1))';
         ys28 = (sum(Y1(1:28,:)))';
         ys=ys28;
         ys(L29)=ys29(L29);
      else
         ys  = (sum(Y1))';
      end
      
   else; % for variables other than pcp, snowfall, #heating deg days and # cooling dd
      % compute mean over days rather than sum
      if n==2; % february special: mean over 29 days inleap yr, 28 in other years
         ys29 = (mean(Y1))';
         ys28 = (mean(Y1(1:28,:)))';
         ys=ys28;
         ys(L29)=ys29(L29);
      else
         ys  = (mean(Y1))';
      end
   end
   Y(:,n+1) = ys;
   
   
end
fclose(fid1);







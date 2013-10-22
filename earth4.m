function earth4
% earth3:  Compute adjusted average minimum daily temperatures for San Pedro stations
% CALL: earth3;
%
% Meko 12-20-97
%
% San Pedro Project.  We needed to produce climatic 'surfaces' of variables,
% including annual average minimum daily temperature.  Susan Skirvin (SS) prepared input file
% of annual minimum ave daily tmp for stations.  I wrote earth4.m to compute the 
% long-term station means, and to adjust the means as needed to be consistent
% with the temporal anomalies in minimum temperature at Tombstone
%
%******************* IN  ********************************
%
% No input arguments
%
% User is prompted for name of file holding the input station data
% produced by SS
%
% User also prompted with chance to change the reference period
%
%********************* OUT **********************************
%
% No output arguments
%
% One ascii file is produced.  The file is hard coded with name adjmint.dat.
% This file has a row of info for each station, and 5 cols:
%
% 1 state -- AZ or NM
% 2 station code
% 3 unadjusted min -- deg F
% 4 adjusted min
% 5 sample size for the min
%
%
%******************  NOTES ***********************************
%
% I decided to use the global HCN urban heat island adjusted Tombstone 
% ave min monthly 
%
% record as a basis for adjustment.  The reference period is 1898-1996.  
% The global HCN for Tombstone does have some missing monthly data in this
% interval. I handle this by filling in the missing monthly values with
% the long-term monthly means of ave daily min tmp computed on whatever
% years in 1898-1996 have data
%
% Related programs are earth1.m, earth2.m and earth3.m.  
% Earth1.m and earth2.m operate on monthly
% data and are vastly different from earth3.m and earth4.m in structure
%
% Earth4.m differs from earth3.m in that I must first estimate some missing
% data in the master series. Another difference is that earth4 uses the global
% Tombstone hcn while earth3 used the US hcn for the master series.  The global
% hcn units are deg C while the US hcn has units deg F
%
% The temperature is adjusted so that the station means are 
% representative of reference-period means.  The adjustment is based
% on the departure from the reference mean for the sample years at
% the master station (Tombstone HCN urban heat island adjusted)
%
%**************************** STEPS **************************
%
% Compute the time series of reference-period-term means for the master series,
%   including the preliminary filling in of missing data
% Convert master reference series from deg C to deg F
%
% Initial read of station data file, assumed to have  4 pieces of data separated by spaces
%  1-state (AZ or NM)
%  2-station number (e.g., 287)
%  3-year
%  4-data value
% Count the number data lines and number of changes in station code.  The 
%   number of changes in station code gives upper limit on number of stations
%   in the filef stations, and number of data lines for each
% Size vectors to store state, station code, year and value
% Second read of data file to store the data
% Loop over stations
%		Compute and store sample mean and adjusted mean
% end loop
% Write output ascii file
% 
%
%
%************************** START CODE ***************************

% Prompt to allow user to change reference period
prompt={'Enter the first year:','Enter the last year:'};
def={'1893','1996'};
tit='Reference Period for Adjustment';
lineNo=1;
answer=inputdlg(prompt,tit,lineNo,def);
yrgo = str2num(answer{1});  yrsp = str2num(answer{2});

% Get the file with the  data for the master series -- global HCN urban ha series
% for Tombstone
pf1='c:\projs\ac4\tmp\susan\e360.mat';
eval(['load ' pf1]); % monthly data is in X; covers 1893-1996;
Z=X; % want data for master in Z

%  Some monthly data are missing for year 1893-1996.  Replace with long-term monthly means
Ztemp = Z(:,2:13);
Ltemp = isnan(Ztemp);
tempmn = nanmean(X(:,2:13));
tempmn = repmat(tempmn,size(Z,1),1);
Ztemp(Ltemp) = tempmn(Ltemp);
Ztemp = (9/5) * Ztemp +32;
Z(:,2:13)=Ztemp;
clear Ltemp Ztemp tempmn;


% Earliest year of any of Susan's 
% Earthinfo series is 1898. I use 1898-1996 as reference period


% Compute annual means of min daily temp from ave monthly means of daily min
Ztemp = Z(:,2:13);
ztemp = (mean(Ztemp'))';
yr = Z(:,1);

% Cull  reference period data
Lref = yr>=yrgo & yr<=yrsp;
Z = ztemp(Lref);
yr = yr(Lref);

% Compute reference-period mean for master series
zm = mean (Z);


% Load Susan's input file with data for all stations
load c:\projs\ac4\tmp\susan\tminday.dat;
W = tminday;
% 4 cols:
%  1- state code (1==AZ, 2==NM)
%  2- station numeric code
%  3- year
%  4- data value (deg F)


%*********  Initial read to count number of files and store names

% Initialize a blank string matrix to hold station codes;
% allow for maximum of 500 stations
b=blanks(4);
B=repmat(b,500,1);

% Initialize years matrix
a=NaN;
N=repmat(NaN,500,1); % will hold number of years

% Size the input data matrix
[mW,nW]=size(W);


% Loop over lines in W
idold = 99999; 
ns= 0 ; % counter for number of stations
ny=0; % counter for number of years of data for a station

for n = 1:mW;
   w = W(n,:);
   if w(2)~=idold; % new station
      ns=ns+1;
      N(ns)=1;
      idold = w(2); % update current id
   else; % another line for this station
      N(ns)=N(ns)+1;
   end
end

      


% Now ns should be number of pseudo-stations and N should hold number
% of years of data for each

% Truncate unused part of N
N=N(1:ns,:);

nlines = sum(N); % total number of lines, not including header

% Initialize vectors to hold data
statenm = repmat(blanks(2),nlines,1);
stncode = repmat(blanks(4),nlines,1);
yrvect = repmat(NaN,nlines,1);
xdata  = repmat(NaN,nlines,1);


%--------------   STORE THE VECTORS
for n = 1:nlines;
   w = W(n,:);
   if w(1)==1;
      statenm(n,:)='AZ';
   elseif w(1)==2;
      statenm(n,:)='NM';
   else
      error('State code should be 1 or 2');
   end
   stncd = int2str(w(2));
   nchar = length(stncd);
   stncode(n,1:nchar) = stncd;
   stncode(n,:) = strjust(stncode(n,:),'right');
end
yrvect=W(:,3);
xdata = W(:,4);



%----------- BUILD LIST OF STATIONS, WITH NO DUPLICATES

state1  = [];
stn1 = [];

ncount=0; % counter for new stations found
for n=1:nlines;
   c1 = statenm(n,:);
   c2 =stncode(n,:);
   if n==1;
      state1=c1;
      stn1=c2;
      ncount=1;
   else
      nold = size(state1,1);
      statec = repmat(c1,nold,1);
      stnc = repmat(c2,nold,1);
      L1 = statec==state1;
      L2 = stnc == stn1;
      
      L3 =    any (all(L1') & all(L2'));
      if ~L3
         ncount=ncount+1;
         state1=[state1; c1];
         stn1=[stn1; c2];
         [state1 stn1];
      end
      
   end
end

nstns = size(state1,1);  % number of uniquely numbered stations

% Initialize matrix to store start, end years
M1=repmat(NaN,nstns,1);  % to hold unadjusted means
M2=repmat(NaN,nstns,1); % to hold adjusted means
S=repmat(NaN,nstns,1);  % to hold sample size for mean


%**************** PART 2: 

% Know this:  have unique station ids in state1 and stn1
%  have station ids for all lines of input file in statenm, stncode
%  have year vector and data in yrvect, xdata

disp('LOOPING OVER STATIONS')

S2=[statenm stncode];

% Loop over stations
for n = 1:nstns;
   
   
   % Find rows of data for this station
   state = state1(n,:);
   stn  = stn1(n,:);
   S1 = [state stn];
   
   %Debugging
   %if strcmp(stn,' 415');
    %  disp('here');
   %end
   
   S1 = repmat(S1,nlines,1);
   L1 = (all((S1 == S2)'))';
   
   N(n)=sum(L1); % Number of years of data at this station
   
   Y = [yrvect(L1) xdata(L1)]; % year and data value
   
   
	% compute the unadjusted mean and adjusted mean
	y = Y(:,2); % pull col of data
	yrstn = Y(:,1);  % vector of years for the station data
	L3=yrstn>=yrgo & yrstn<=yrsp;  % logical pointer to years in master coverage
	v = y(L3);   % col vect of station data in reference period
	t=yrstn(L3); % col vector of corresp years station has data in ref  period
	mn1=mean(v);  % unadjusted mean
	ss1=sum(L3); % sample size for that mean
		
	% Get the corresp years of the master series and compute 
	% its subperiod mean
	iz = t-yrgo+1; % pointer to rows of Z
	z =Z(iz); % cv of data for master series
	zms = mean(z); % sub-period mean
     
	
	%*********** Compute Adjusted Mean  -- TMP ****************
	r= zms-zm; % difference of subperiod mean temperature
				% and full period mean temperature for master
	% Compute adjusted mean for the sample
	mn2=mn1-r;
      
   %************************************************************
	% Handle the special case in which the current station has
	% at least one year of data before start of the master series
	% or after end of master series
   L2 = yrstn<yrgo | yrstn>yrsp;
   if any(L2);  % at least one year's data outside the reference period
      fclose all;
      disp(['Station ' B(n,:) ]);
      error('Above series has data outside reference period');
      
      % For earth3, I do not us the following dcode in this if segment
		voutside = y(L2);   % col vect of data outside ref pd
		toutside=yrstn(L2);  % col vector of corresp years
		mnout=mean(voutside);  % unadjusted mean for data outside ref pd
		sout=sum(L2); % sample size for outside-ref period
		% Compute weighted mean
		mn1=((sout*mnout)+(mn1*ss1))/(sout+ss1);
		mn2=((sout*mnout)+(mn2*ss1))/(sout+ss1);
		ss1=ss1+sout; % will store total sample size
	else;
	end

   
   % Store the Unadjusted Mean, Adjusted Mean, and Sample size
	M1(n)=mn1; % unadjusted
	M2(n)=mn2; % adjusted
	S(n)=ss1; % sample size

   
end

N = N(1:nstns);


%********************* PART 3: MAKE ASCII FILES OF OUTPUT

disp('MAKING ASCII OUTPUT FILES');

file4='adjmint.dat';

fid4=fopen(file4,'w');
fmt1 = '%10.5f';
fmt2 = '%5.0f';

header='station      Unadj       Adj    N(yr)';
fprintf(fid4,'%s\n',header');

for n = 1:nstns; % loop over stations
   str1 = sprintf('%2s %6s',state1(n,:),stn1(n,:));
	str2=sprintf(fmt1,M1(n)); % unadj means
	str3=sprintf(fmt1,M2(n)); % adj means
	str4=sprintf(fmt2,S(n)); % sample sizes

	fprintf(fid4,'%s %s %s %s\n',str1,str2, str3, str4);
end


fclose all;


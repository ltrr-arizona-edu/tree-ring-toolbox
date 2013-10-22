function earth3
% earth3:  Compute adjusted annual mean temperatures for San Pedro stations
% CALL: earth3;
%
% Meko 12-20-97
%
% San Pedro Project.  We needed to produce climatic 'surfaces' of variables,
% including annual mean temperature.  Susan Skirvin (SS) prepared input file
% of annual mean time series for stations.  I wrote earth3.m to compute the 
% long-term station means, and to adjust the means as needed to be consistent
% with the temporal anomalies in annual mean temperature at Tombstone
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
% One ascii file is produced.  The file is hard coded with name adjannt.dat.
% This file has a row of info for each station, and 5 cols:
%
% 1 state -- AZ or NM
% 2 station code
% 3 unadjusted mean -- deg F
% 4 adjusted mean
% 5 sample size for the mean
%
%
%******************  NOTES ***********************************
%
% I decided to use the HCN urban heat island adjusted Tombstone monthly 
% record as a basis for adjustment.  The reference period is preferably
% 1893-1996, which has no missing monthly values
%
% Related programs are earth1.m and earth2.m, which operate on monthly
% data and are vastly different from earth3.m in structure
%
% The temperature is adjusted so that the station means are 
% representative of reference-period means.  The adjustment is based
% on the departure from the reference mean for the sample years at
% the master station (Tombstone HCN urban heat island adjusted)
%
%**************************** STEPS **************************
%
% Compute the time series of reference-period-term means for the master series
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
%*******************  NOTES ***********************************
%
% Special case -CODE COMMENTED OUT IN EARTH3.
% say a station has years before or after the perio
% covered by the master.  Then I compute a weighted mean, weighted
% by subsample size.  For the years in the master coverage I adjust.
% For the earlier years, I do not adjust, but compute a straight mean.
% Then I use  m = (N1*m1 + N2*m2)/(N1+N2) for weighted mean. For 
% example, say master covers 1898-96, and a station has data for
% 1892-1970.  Compute m1 on 1898-1970, adjusted. Compute m2 on 
% 1892-1897.  compute weighted mean as
%
% m = (73*m1+6*m2)/(73+6)
%
% Note that I have commented out this code for earth3.m.  Earth3.m
% simply returns and error message if any station has data outside
% the reference period.  No problem with annual T for this data set
% because no stations have such data (if ref period 1893-1996)
%
%************************** START CODE ***************************

% Prompt to allow user to change reference period
prompt={'Enter the first year:','Enter the last year:'};
def={1893,1996};
tit='Reference Period for Adjustment';
lineNo=1;
answer=inputdlg(prompt,tit,lineNo,def);
yrgo = str2num(answer{1});  yrsp = str2num(answer{2});

% Get the file with the  data for the master series -- HCN urban ha series
%[file1,path1]=uigetfile('mast*.dat','File with master series data');
pf1='c:\projs\ac4\tmp\monhcnu\tombt.mat';
eval(['load ' pf1]); % monthly data is in Z; covers 1889-1996;
%  No monthly data are missing for year 1893-1996;
% Earthinfo series starts in 1898. I use 1898-1996 as reference period


% Compute annual mean temp from monthly data
Ztemp = Z(:,2:13);
ztemp = (mean(Ztemp'))';
yr = Z(:,1);

% Cull  reference period data
Lref = yr>=yrgo & yr<=yrsp;
Z = ztemp(Lref);
yr = yr(Lref);

% Compute reference-period mean for master series
zm = mean (Z);


% Get the file with the monthly data for the other stations
[file2,path2]=uigetfile('avg_all.dat','File with earthinfo data');
pf2=[path2 file2];
fid2=fopen(pf2,'r');


%*********  Initial read to count number of files and store names

% Initialize a blank string matrix to hold station codes;
% allow for maximum of 500 stations
b=blanks(4);
B=b(ones(500,1),:);

% Initialize years matrix
a=NaN;
N=a(ones(500,1),:); % will hold number of years


% Hard code the columns that expect the 4 pieces of info to be in 
c1go=1; c1sp=2; % state code
c2go=17; c2sp=20; % id
c3go=28; c3sp=31; % year
c4go=44; c4sp=51; % data value


clc
disp('INITIAL READ TO COUNT LINES');
% Loop over lines in the file
c=fgetl(fid2); % skip header line
idold = '999'; % initialize old id
ns= 0 ; % counter for number of stations
ny=0; % counter for number of years of data for a station
k1=1; % control for while
while k1==1
	c=fgetl(fid2); % get line of file
	if feof(fid2) & length(c)<2;
		N(ns)=ny; % store number of years of data
		k1=0;
	else
		id =c(17:20);
		if ~strcmp(id,idold); % new station
			if ns~=0; % Not first station
				N(ns)=ny; % store number of years of data in prev stn 
			end
			ns=ns+1; % increment station counter
			B(ns,1:4)=id;	% store 3-char station code	
			ny=1; % re-set counter for years for this station
         %disp(['Station count  = ',int2str(ns)]);
         idold=id;
      else
         idold=id;
			ny=ny+1; % increment year counter
		end
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


%-------------- SECOND READ OF FILE2:  STORE THE VECTORS
disp('SECOND READ, TO STORE STRING AND NUMBERIC INFO IN MATRICES');

frewind(fid2);
fgetl(fid2); % skip header line
for n = 1:nlines;
   c=fgetl(fid2);
   statenm(n,:) = c(c1go:c1sp);
   stncode(n,:)=c(c2go:c2sp);
   yrvect(n)= str2num(c(c3go:c3sp));
   xdata(n)=str2num(c(c4go:c4sp));
end


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


% Rewind the input file
fclose(fid2);


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

file4='adjannt.dat';

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
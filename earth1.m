function earth1
%
% Compute time-adjusted monthly means from earth-info generated
% monthly pcp data.
%
% Meko 8-18-96
%  rev 9-10-96 to handle case of a station/month having some data
%		before start of master and all that early data missing values
%
%******************* IN ARGS ********************************
%
% None.  Two files must be available for user to click. One
% has the Earth-Info file of multi-station monthly data. The other
% has the regional-mean monthly data 
%
% This function makes three output files with hard-coded names.
%
% 1 mnorig.dat   original means -- before adjustment
% 2 mnadj.dat	adjusted means
% 3 nsize.dat  sample size for means
%
% The output files of means are ascii.  Each line has
% the 5-character station code followed by the 12 values
% in inches.
%
%**************************** STEPS **************************
%
% Compute the long-term monthly means for the master regional series
% Read the data file, assumed to have *ccccc starting each new stn
%	Count the number of stations, and number of data lines for each
% Rewind the data file
% Loop over stations
%		Initialize a storage matrix for the monthly series
%		Read file line by line
%			Convert blank to NaN
%			Compute row-index for storage from year
%			Store the year's data as numeric in matrix
%			Reach end year for this station
%		Operate on the numeric matrix
%			Loop over months 1-12
%				Pull subset of valid rows (not NaN)
%				Compute unadjusted means from this vector
%				From the years for the vector, compute a row index
%					into master series
%				Pull the corresp subset of rows of master series
%				Compute the mean for the master series, and
%					the ratio of that mean to the long-term
%					master series mean
%				Adjust the mean for the current series		
%		Store the unadusted and adjusted means for this station
%
%
%
%*******************  NOTES ***********************************
%
% Special case -- say a station has years before or after the perio
%	covered by the master.  Then I compute a weighted mean, weighted
% by subsample size.  For the years in the master coverage I adjust.
% For the earlier years, I do not adjust, but compute a straight mean.
% Then I use  m = (N1*m1 + N2*m2)/(N1+N2) for weighted mean. For 
% example, say master covers 1904-95, and a station has data for
% 1895-1970.  Compute m1 on 1904-1970, adjusted. Compute m2 on 
% 1895-1903.  compute weighted mean as
%
% m = (67*m1+9*m2)/(67+9)
%
%
% Special case:  a stations has a short series.  So short that 
% maybe one of the means for the master series for the sub-period
% is zero.  This would lead to a zero means ratio, and infinity
% when dividing the sample-station mean by the ratio.  In this case,
% just accept the sample-station mean.

%************************** START CODE ***************************


% Get the file with the monthly data for the master series
[file1,path1]=uigetfile('mast*.dat','File with master series monthly');
pf1=[path1 file1];
eval(['load ' pf1]);

% Put the master series into  Z (12 cols), and year into yr (1 col)
% Compute the long-term means for the master series
f1=strtok(file1,'.');
eval(['Z = ' f1 ';'])
yr = Z(:,1);
yrgo= yr(1); % start year of mast er series
Z=Z(:,2:13);
zm=mean(Z); % long-term means (row vector, 12 values) of master
if any(isnan(zm)),
	error('Master series had at least one NaN')
end

% Get the file with the monthly data for the other stations
[file2,path2]=uigetfile('reg*.txt','File with earthinfo data');
pf2=[path2 file2];
fid2=fopen(pf2,'r');


%*********  Initial read to count number of files and store names

% Initialize a blank string matrix to hold station codes;
% allow for maximum of 500 stations
b=blanks(5);
B=b(ones(500,1),:);

% Initialize years matrix
a=NaN;
N=a(ones(500,1),:); % will hold number of years

% Loop over lines in the file
ns= 0 ; % counter for number of stations
ny=0; % counter for number of years of data for a station
k1=1; % control for while
while k1==1
	c=fgetl(fid2); % get line of file
	if feof(fid2)
		N(ns)=ny; % store number of years of data
		k1=0;
	else
		if c(1)=='*'; % new station
			if ns~=0;
				N(ns)=ny; % store number of years of data in prev stn 
			end
			ns=ns+1; % increment station counter
			B(ns,1:5)=c(2:6);	% store 5-char station code	
			ny=0; % zero the counter for years for this station
			disp(['Station count  = ',int2str(ns)]);
		else
			ny=ny+1; % increment year counter
		end
	end
end


% Now ns should be number of stations and N should hold number
% of years of data for each

% Truncate unused part of N
N=N(1:ns,:);

% Initialize mtarix to store start, end years
YRS = a(ones(ns,1),ones(2,1)); % will hold start, end years
M1=a(ones(ns,1),ones(12,1)); % to hold unadjusted means
M2=a(ones(ns,1),ones(12,1)); % to hold adjusted means
S=a(ones(ns,1),ones(12,1)); % to hold sample size for mean


% Rewind the input file
frewind(fid2);


%**************** PART 2: 


clc
disp('STARTING PART 2: DATA MANIPULATION')

d=a(:,ones(12,1)); % will store monthly data values for a year

for n=1:ns; % loop over stations
	disp(['Part 2; working on station: ' int2str(n)]);
	%	Initialize a storage matrix for the monthly series
	nc = N(n); % station should have this many years of data
	Y = a(ones(nc,1),ones(13,1)); % will hold the year and monthly data
	c = fgetl(fid2); % get a line
	if c(1)~='*',
		error('First line of block for stn must start with *')
	end
	% Loop over the lines of data 
	for j =1:nc;
		c=fgetl(fid2);
		if c(1)=='*';
			error('Expected data line, got *')
		end
		year=str2num(c(9:12)); % year
		for i=1:12; % loop over the monthly data vlues
			ngo= 13+(i-1)*6; % start index for field
			nsp=ngo+5; % end index
			cfield= c(ngo:nsp);
			if all(isspace(cfield)); % if the file is blank
				d(i)=NaN;
			else
				d(i)=str2num(cfield);
			end
		end
		% Store this year's year and 12 values
		Y(j,:)=[year d];
	end

	% Have matrix Y with year column and 12 monthly values, some
	% maybe NaN.

	% Troubleshooting segment
	if n==31; % if Douglas
		disp(['n = ' int2str(n)]);
	end

	
	% Loop over the 12 months of the year for this station to
	% compute the unadjusted mean and adjusted mean
	for i=1:12; % loop over months
		y = Y(:,i+1); % pull col of monthly data
		L1=~isnan(y); % pointer to valid rows
		L2=Y(:,1)>=yrgo;  % years in master coverage
		L3=L1 & L2;  % valid data and in master series period
		v = y(L3);   % col vect of valid data; also must be in period
			% covered by master series
		t=Y(L3,1); % col vector of corresp years
		mn1=mean(v);  % unadjusted mean
		ss1=sum(L3); % sample size for that mean
		
		% Get the corresp years of the master series and compute 
		% its subperiod mean
		iz = t-yrgo+1; % pointer to rows of Z
		z = Z(iz,i); % cv of data for master series
		zms = mean(z); % sub-period mean
		% Handle the special case of weird sub-period with sub-period
		% monthly mean zero
		if zms>0; % if sub-period mean for master greater than zero
			r = zms/zm(i); % ratio of sub-period mean to full mean for
				% master series
		else
			r=1; % to avoid dividing by zero later
		end
		
		% Adjust the mean for the sample
		mn2=mn1/r;  % adjusted mean

		% Handle the special case in which the current station has
		% at least one month of valid data before start 
		% of the master series
		if ~all(L2);  % at least one year's data before start of master
			L4=L1 & ~L2; % point to valid years before start of master
			if sum(L4)==0; % all the earlier values are NaN; no action needed
				disp(['Station ' int2str(n)]);
				disp('Rare case: data before start of master all NaN');
			else; % some valid data before start of master;  must compute
				% its mean
				vearly = y(L4);   % col vect of valid early data
				tearly=Y(L4,1); % col vector of corresp years
				mnearly=mean(vearly);  % unadjusted mean for ealy years
				searly=sum(L4); % sample size for that early mean
				% Compute weighted mean
				mn1=((searly*mnearly)+(mn1*ss1))/(searly+ss1);
				mn2=((searly*mnearly)+(mn2*ss1))/(searly+ss1);
				ss1=ss1+searly; % will store total sample size, early and late
			end
		end

		% Store the Unadjusted Mean, Adjusted Mean, and Sample size
		M1(n,i)=mn1; % unadjusted
		M2(n,i)=mn2; % adjusted
		S(n,i)=ss1; % sample size
	end; % of loop over months

end; % of loop over stations

fclose (fid2);

%********************* PART 3: MAKE ASCII FILES OF OUTPUT
clc
disp('STARTING PART 3: MAKING ASCII OUTPUT FILES');

file4='mnorig.dat';
file5='mnadj.dat';
file6='nsize.dat';

fid4=fopen(file4,'w');
fid5=fopen(file5,'w');
fid6=fopen(file6,'w');

fmt1='%7.3f'; % format for means
fmt2='%5.0f'; % format for sample size

for n = 1:ns; % loop over stations
	str1=sprintf('%s',B(n,:)); % code
	str2=sprintf(fmt1,M1(n,:)); % unadj means
	str3=sprintf(fmt1,M2(n,:)); % adj means
	str4=sprintf(fmt2,S(n,:)); % sample sizes

	fprintf(fid4,'%s %s\n',str1,str2);
	fprintf(fid5,'%s %s\n',str1,str3);
	fprintf(fid6,'%s %s\n',str1,str4);
end



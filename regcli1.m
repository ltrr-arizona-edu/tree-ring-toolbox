function regcli1
%
% Regional PPT or TMP series from station data by standardized normal method.
% Produces one regional series only (not series at gridpoints)
% Weights standardized anomalies at stations evenly (i.e., no 
%		inverse-distance or other type of weighting scheme used
% 
% Meko 1/17/97
%
%***************************** INPUT FILES *********************************
%
% 1. fls*.txt -- edix-produced file with list of prefixes of filenames
%		of .mat files of the raw monthly climate data.  First line
%		should give path. For example:
%
%			d:\sprb\ppt\rawmon\
%			tucs 1
%			bisb 2
%			doug 3
%			...etc
%
%		The number of lines, or filenames, in this file is ns1. The
% 		number to the right of the file name is ignored by the program
%
% 2. misc*.mat -- miscellaneous program control
%		im (m1 x 1)r pointer to rows of file1 telling which stations to
%			use for the master series
%		ir (m2 x 1)r pointer to rows of file1 telling which to use for the
%			regional series
%		pdref (1 x 2)i reference period (start, end years) for adjusting
%			means and standard deviations
%		pdstore (1 x 2)i  start, end years for storage matrix that must
%			cover at least the period with any valid data at any station.
%		kmode (1 x 1)i mode of run.  kmode==1 means diagnostic. kmode==2
%			means skip diagnostic output. Diagnostic run gives info useful
%			for debugging.  (a) sceen output for stations showing segment
%			of monthly data, means and standard deviations, and standardized
%			normal data.
%		kopts (1 x 2)i options
%			kopts(1):  	==1 --> do not adjust std deviations, just means
%					==2 --> adjust means and std deviations
%			kopts(2):	==1 --> subst. zero for negative monthly regl values
%					==2 --> accept negative monthly regional values
%
%					(note: typically will set ==1 for ppt, ==2 for temper.)
%**************************** OUTPUT FILES ***********************************
%
% 1. umean*.txt ascii file of unadjusted monthly means for each stations.
%			Each row contains station identifier from file name and 12 values
%			Last line is mean over stations in region
% 2. amean*.txt similarly, adjusted means
% 3. ustdev*.txt ... unadjusted standard devs
% 4. astdev*.txt ... adjusted standard devs
% 5. usize*.txt ... sample size (# years in reference period) for unadjusted
%		means and st devs
% 6. uregsn*.dat regional standardized normal time
%		series matrix, version based on no adjustment of means and st devs
% 7  aregsn*.dat like uregsn*.dat, but basedon adjusted means and st devs
% 8  nstns*.dat time series matrix of number of stations in regional average
% 9  ureg*.dat regional climate series, in original data units, based on
%		analysis without adjusting means and st devs
% 10  ureg*.mat ... like ureg*.dat, but a .mat file
% 11 areg*.dat regional climate series,..., analysis with adjustment of
%			means and std devs.
% 12 areg*.mat ... like areg*.dat, but a .mat file
% 13 ratmn*.dat monthly ratios used to adjust station means
% 14 ratsd*.dat monthly ratios used to adjust station std devs
%
%******************************** NOTES *************************************
%
% Reference: Jones & Hulme (1996).  The reference does not cover the case
% of adjusting the means and standard deviations of station data to
% account for anomalous climatological reference periods.
%
% Method.   Monthly station time series are first converted to standardized
% normal series using the mean and standard deviation for that station/month.
% The standardized normal series are then averaged over stations to produce
% a regional standardized normal series.  The regional series is converted
% from standardized units to original units (e.g., mm of PPT) using the
% a regional monthly mean and regional monthly standard deviation. 
%
% The monthly means and standard deviations are optionally adjusted to a
% common reference period.  This might be desireable when the period of 
% available data for one or more stations has anomalous mean or standard 
% deviation relative to the long-term mean and standard deviation.  The 
% adjustment is done by ratios using a "master" subregional series built
% from stations with long records.  Consider a "biased" series and a 
% master series.  The means and standard deviations for the master series 
% are computed for the specified long-term reference period and for the
% subset of years with data at the biased station.  The ratios of the
% long-term and subperiod means at the master series is computed for 
% Jan, Feb, ..., Dec.  Ratios are similarly computed for the standard deviation.
% The means and standard deviations for the biased series are then
% adjusted by multiplying them by the computed ratios.
%
%******************************************************************************

% Allocate
blnks=blanks(8);
B1=blnks(ones(200,1),:); % Initialize to allow for 200 filenames

%********************************************************************
% Get list of input .mat files and store in flnames.  Also compute number
% of .mat files in the list.  First line of the list is the path to the files
[file1,path1]=uigetfile('fls*.txt','List of input .mat files');
pf1=[path1 file1];
fid1=fopen(pf1,'r');
path3=strtok(fgetl(fid1));
k1=1;
ns1=0; % counter for number of .mat files in list
while k1;
	c1=fgetl(fid1);
	if feof(fid1);
		k1=0;
	else
		if isspace(c1(1));
			error('Filenames in list must start in col 1');
		end
		ns1=ns1+1;
		c2=strtok(c1);
		nc2=length(c2);
		B1(ns1,1:nc2)=c2;
	end
end
B1=B1(1:ns1,:);

% B1 now contains file name prefixes
% ns1 is the number of files
% path3 is the path to the .mat input station data files


%*********** Read and QC  miscellaneous control information

[file2,path2]=uigetfile('misc*.mat','Miscellaneous program control');
pf2=[path2 file2];
eval(['load ' pf2]);

% Check that all variables wanted from the misc file exist
Ltemp1=[exist('im') exist('ir') exist('pdref')];
Ltemp2=[exist('kmode') exist('pdstore') exist('kopts')];
Ltemp=[Ltemp1 Ltemp2];
Ltemp=all(Ltemp==1);
if ~Ltemp
	clc
	disp('The following variables should be in the misc.mat file:');
	disp('im,ir,pdref,pdstore,kmode,kopts');
	error('All the above were not in the file');
end

% Check that values of im, ir, pdref, kmode acceptable
if max(im)>ns1 | max(ir)>ns1;
	error('Some entry in im or ir too large for number of listed stns');
end
if pdref(1)>pdref(2); 
	error('Reference period ends before it starts');
end
if pdstore(1)>pdstore(2)
	error('specified values for pdstore impossible');
end
if kmode<1 | kmode>2,
	error('Invalid setting for kmode')
end

if kopts(1)<1 | kopts(1)>2;
	error('Illegal setting of kopts(1)')
end
if kopts(2)<1 | kopts(2)>2;
	error('Illegal setting of kopts(2)')
end


% Calc number of stations in regional series and in master
ns2=length(ir);
ns3=length(im);

%Display on the screen the files specified for master and regional makeup
B2=B1(im,:);
B3=B1(ir,:);
clc;
disp('Series specified for master')
stemp=sprintf('%\n');
disp(B2);
disp(stemp)
disp('Series specified for region')
disp(B3);
disp(stemp)
disp('Press any key to continue')
pause

%****************************************************************
% Compute the master series for the reference period, and the 
% reference-period monthly means and standard deviations for the
% master series. 
% Note that already know that
% ns3 is number of stations in the master series
% B2  holds file names of the master series stations
disp('COMPUTING MASTER SERIES AND RELATED QUANTITIES')

% Size storage matrices
a=NaN;
yr1=(pdref(1):pdref(2))';
nref=length(yr1); % # yrs in ref pd
ncols=ns3*12; % number of columns needed
M=a(ones(nref,1),ones(ncols,1)); % store ref pd raw data, years by months
LL=a(ones(nref,1),ones(ncols,1)); % store ref pd missing value pointer
V=a(ones(nref,1),ones(ncols,1)); % store ref period station stdzd anomalies
VM=a(ones(nref,1),ones(12,1)); % mean stdzd anomaly series for the master
YM=a(ones(nref,1),ones(12,1)); % VM converted to original climate units
Mrefm=a(ones(ns3,1),ones(12,1)); % ref period long-term monthly means
Mrefs=a(ones(ns3,1),ones(12,1)); % ref period long term std devs
gm=a(:,ones(12,1)); % global mean monthly ppt (average over stations)
gs=a(:,ones(12,1)); % global mean std dev

% Get the monthly data for the master stations and store in matrix M
for n = 1:ns3;
	fln=strtok(B2(n,:));
	eval(['load ' path3 fln]);
	yr = Z(:,1);
	L1=yr>=pdref(1) & yr<=pdref(2);
	Z=Z(L1,2:13);  % pull subset of rows
	if size(Z,1)~= nref
		disp(fln)
		error('Above station for master does not cover reference period');
	end
	% Compute the long-term monthly means and std deviations
	L1= isnan(Z); % Pointer to NaNs in Z
	if any(any(L1)); 
		[zmn,nsize]=meannan(Z);
		[zst,nsize]=stdnan(Z);
	else; % no missing data in any months
		zmn=mean(Z);
		zst=std(Z);
	end

	% Compute standardized anomalies
	ZMN=zmn(ones(nref,1),:);
	ZST=zst(ones(nref,1),:);
	V1=(Z-ZMN) ./ ZST;

	% Compute storage columns
	jgo = (n-1)*12 +1;
	jsp = jgo+11;
	jcols=[jgo:jsp];

	% Store results
	M(:,jcols)=Z; % monthly data
	V(:,jcols)=V1; % standardized anomalies
	Mrefm(n,:)=zmn; % ref-pd monthly means
	Mrefs(n,:)=zst; % ref-pd monthly std devs
end

% Compute global mean and std dev
gm=mean(Mrefm);
gs=mean(Mrefs);

% Compute ref-pd time series of  standardized anomalies averaged over master
% series
S1=zeros(nref,12); % allocate for a sum
for m=1:nref; % loop over ref period years
	for k = 1:12; % over months
		j1=(1:12);
		J1=j1(ones(ns3,1),:);
		j2=(1:ns3)';
		j3=12*(j2-1);
		J3=j3(:,ones(12,1));
		J4=J1+J3;
		jcols=(J4(:,k))';
		v=V(m,jcols);
		v(isnan(v))=[];
		vm=mean(v); % mean standardized anomaly, this month/year
		VM(m,k)=vm; % tsm of mean standardized anomalies
	end
end

% Convert the master series of mean standardized anomalies to climate
% series in original units
GM=gm(ones(nref,1),:); % expand rv of global means to a matrix
GS=gs(ones(nref,1),:); % expand.. std devs
YM=GM + GS .* VM;

% Store the reference-period standard deviations for the master series. 
% These 12 values are not the same as gs because the mean of the standard
% deviations of several series is not the same as the standard deviation
% of the mean series.
gsnew=std(YM);


%
%  ALL DONE WITH GETTING "MASTER-SERIES" REFERENCE PERIOD MEANS AND STD DEVS
%
% YM tsm of regional climate variable for the master series, in original units
%		for the reference period
% gm rv of master-series monthly means for reference period
% gsnew rv of master-series monthly std devs for ref pd
%
%
%**************************************************************************
%
disp('READING MONTHLY DATA FILES AND STORING DATA')
% Store the input monthly ppt series for all specified climate stations for
% the region in a single matrix. Note that:
%
% ir -- pointer to rows of B3 telling stations to form region 
% B3 -- filename matrix with those stations' file names
% ns2 -- number of stations for regions (equals length(ir))
% pdstore -- start, end years for multistation multimonthly storage matrix
% pdref -- start, end years of reference period
% yr1 -- cv of years for reference period
% nref -- number of years in reference period

% Initialize matrices 
yr2=(pdstore(1):pdstore(2))'; % years vector for storage matrix
n2 = length(yr2); % number of rows (years) in storage matrix
ncols = ns2*12; % % number of cols in storage matrix
D=a(ones(n2,1),ones(ncols,1));  % monthly raw data
E=a(ones(n2,1),ones(ncols,1));  % monthly stdzd anomalies, based on
%		unadjusted means and std devs
F=a(ones(n2,1),ones(ncols,1));  % monthly stdzd anomalies, based on
%		adjusted means and std devs
T=a(ones(ns2,1),ones(2,1));  % start, end years of each station's record
I=a(ones(ns2,1),ones(2,1));  % corresp row pointer to D,E


% What cols in D and E will the station monthly data go into?
j1=(1:12:(ncols-11)); % start col for stn1, 2,etc
j2=j1+11; % end col
J1=[j1' j2'];
% First row of J1 will give first and last storage row for stn 1, second row
%		stn 2 etc

% What cols in D and E will hold all Jan data, Feb data, etc?
j1=(1:12);
j1=j1(ones(ns2,1),:);
j2=(1:ns2)';
j2=12*(j2-1);
j2=j2(:,ones(12,1));
J2=j1+j2;
% First col of J2 will point to cols of D for Jan, secnd for Feb, etc


% Loop over stations going into regional series: get first and last year,
% 
for n=1:ns2;
	clear Z;
	fln=strtok(B3(n,:));
	eval(['load ' path3 fln]);
	if (exist('Z')~=1)
		error(['No Z in: ' fln]);
	end

	% Make year vector; store first year and last year; calc and store
	% row reference info
	yr = Z(:,1);
	nyrs=length(yr);
	if yr(1)<pdstore(1) | yr(nyrs)>pdstore(2)
		disp(fln)
		error('pdstore inconsistent with this series years');
	end
	T(n,:)=[yr(1) yr(nyrs)];
	I(n,:)=[yr(1)-pdstore(1)+1  yr(nyrs)-pdstore(1)+1];


	% Store the monthly data
	D(I(n,1):I(n,2),J1(n,1):J1(n,2)) = Z(:,2:13);
end

%******************************************************************
%
% Pull the monthly data for reference period for each station.
% Make logical pointer matrix to non-nan data in each month
% Store valid sample size for ref period
% Compute the reference-period monthly means for the station
%		also the std devs; store
% Use pointer to loop over months and compute sub-ref-pd monthly
%		means for the master station; same for st devs
% Compute adjustment ratios for means and std devs
% Compute adjusted station means and std devs and store
disp('COMPUTING ADJUSTED MEANS AND STANDARD DEVIATIONS FOR STATIONS')

umean=a(ones(ns2,1),ones(12,1));
ustdev=a(ones(ns2,1),ones(12,1));
usize=a(ones(ns2,1),ones(12,1));
amean=a(ones(ns2,1),ones(12,1));
astdev=a(ones(ns2,1),ones(12,1));
ratmn=a(ones(ns2,1),ones(12,1));
ratsd=a(ones(ns2,1),ones(12,1));

igo=pdref(1)-pdstore(1)+1;
isp=pdref(2)-pdstore(1)+1;
irows=(igo:isp);

for n=1:ns2;
	jcols=J1(n,1):J1(n,2);
	X=	D(irows,jcols); % reference period monthly data
	L1=~isnan(X);
	usize(n,:)=sum(L1);  % store sample size (yrs in ref pd with valid data)
	[mn1,nn]=meannan(X); % unadjusted ref pd means
	[std1,nn]=stdnan(X); % unadjusted std devs
	umean(n,:)=mn1;
	ustdev(n,:)=std1;

	% Master series means and std devs for the sub period of ref period
	L2=~L1; % pointer to NaNs in station series
	sum1=sum(sum(L2)); 
	YYM=YM;
	YYM(L2)=a(ones(sum1,1),:); 
	[mny,nn]=meannan(YYM);
	[stdy,nn]=stdnan(YYM);

	% Ratios for adjusting station means and std devs
	ratmn(n,:)=gm ./ mny;
	if kopts(1)==2; % adjustd std devs
		ratsd(n,:)=gsnew ./ stdy;
	elseif kopts(1)==1; % do not adjust
		ratsd(n,:)=ones(1,12);
	else
		error('Invalid setting for kopts(1)')
	end

	% Adjust means and std devs
	amean(n,:)=ratmn(n,:) .* mn1;
	astdev(n,:)=ratsd(n,:) .* std1;
end




%********************************************************************
%
% Compute and store regional average means and standard deviations
ustats=a(ones(2,1),ones(12,1)); % storage for means in row 1, sdevs in 2
astats=a(ones(2,1),ones(12,1)); % storage for means in row 1, sdevs in 2
% Unadjusted 
ustats(1,:)=mean(umean);
ustats(2,:)=mean(ustdev);
astats(1,:)=mean(amean);
astats(2,:)=mean(astdev);



%*******************************************************************
%
% Fill matrices E and F with standardized departure series for
% stations; 
disp('COMPUTING REGIONAL SERIES')

for n=1:ns2; % Loop over stations in region
	irows=I(n,1):I(n,2); % points to rows in D holding this station's data
	jcols=J1(n,1):J1(n,2); % ...cols in D ....
	X=D(irows,jcols); % Monthly data, entire data record for this station
	[mX,nX]=size(X);

	% Using unadjusted stats
	mns = umean(n,:);
	sds = ustdev(n,:);
	MNS=mns(ones(mX,1),:); % rv to matrix
	SDS=sds(ones(mX,1),:);
	E(irows,jcols)= (X - MNS) ./ SDS;

	% Using adjusted stats
	mns = amean(n,:);
	sds = astdev(n,:);
	MNS=mns(ones(mX,1),:); % rv to matrix
	SDS=sds(ones(mX,1),:);
	F(irows,jcols)= (X - MNS) ./ SDS;
end

%***************************************************************
%
% Compute and store the regional mean standardized departure series,
% unadjusted and adjusted


% Size storage matrices to cover inclusive data period
L1=(any((~isnan(E))'))';
yr3=yr2(L1);
nyrs3=length(yr3);
ZU=a(ones(nyrs3,1),ones(12,1)); % to hold unadjusted regional series
ZA=a(ones(nyrs3,1),ones(12,1)); % to hold adjusted
NU=a(ones(nyrs3,1),ones(12,1)); % to hold unadjusted sample size
	%(number of stations)
NA=a(ones(nyrs3,1),ones(12,1)); % to hold adjusted sample size

i2rows=find(L1); % rows in E to get

% Loop over the 12 months of the year
for k=1:12;
	jcols = J2(:,k); % these columns for this month

	% unadjusted version
	G=E(i2rows,jcols);
	[g1,ng1]=meannan(G');
	ZU(:,k)=g1';
	NU(:,k)=ng1';

	% adjusted version
	G=F(i2rows,jcols);
	[g1,ng1]=meannan(G');
	ZA(:,k)=g1';
	NA(:,k)=ng1';

end

%***********************************************************************
%
% Convert the regional standardized departures to original units of the
% climate variable
% Recall that have unadjusted global means and std devs in ustats, and
% for adjusted series in astats.  
% Recall that have regional standardized departures in ZU, ZA


% Unadjusted first 

% 1x12 vectors of global means and std devs duped to matrices
mu=ustats(1,:);
MU=mu(ones(nyrs3,1),:);
su=ustats(2,:);
SU=su(ones(nyrs3,1),:);

% Stdzd departure to original units
YU=MU + (ZU .* SU);

% Adjusted second 

% 1x12 vectors of global means and std devs duped to matrices
ma=astats(1,:);
MA=ma(ones(nyrs3,1),:);
sa=astats(2,:);
SA=sa(ones(nyrs3,1),:);

% Stdzd departure to original units
YA=MA + (ZA .* SA);


% Optionally substitute zero for negative monthly regional values
if kopts(2)==1; % substitute -- most appropriate for ppt data
	L3=YU<=0;
	numneg=sum(sum(L3));
	if numneg>0;
	 	YU(L3)=zeros(numneg,1);
	end
	L3=YA<=0;
	numneg=sum(sum(L3));
	if numneg>0;
		YA(L3)=zeros(numneg,1);
	end
elseif kopts(2)==2; % do not substitute -- most approp. for tmp data
	% no action
else
	error('Invalid setting for kopts(2)');
end

%******************** BUILD OUTPUT FILES ******************************
%
disp('BUILDING OUTPUT FILES')
 
% Number suffix to be applied to all files
numsuf=input('Number suffix for out files: ');
strsuf=int2str(numsuf);


% The unadjusted reference-period monthly means for stations 
% and the global average over stations

file3=['umean' strsuf '.txt'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt1a='%s ';
fmt1e='%s\n';
fmt1b='%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f';
fmt1c='%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n';
fmt1d=[fmt1b fmt1c];
	
for n=1:ns2;
	b3=B3(n,:);
	x=umean(n,:);
	fprintf(fid3,fmt1a,b3);
	fprintf(fid3,fmt1d,x);
end
fprintf(fid3,fmt1e,' '); % Blank line
fprintf(fid3,fmt1a,'Ave     ');
fprintf(fid3,fmt1d,ustats(1,:));

fclose(fid3);


% The adjusted reference-period monthly means for stations 
% and the global average over stations

file3=['amean' strsuf '.txt'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt1a='%s ';
fmt1e='%s\n';
fmt1b='%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f';
fmt1c='%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n';
fmt1d=[fmt1b fmt1c];
	
for n=1:ns2;
	b3=B3(n,:);
	x=amean(n,:);
	fprintf(fid3,fmt1a,b3);
	fprintf(fid3,fmt1d,x);
end
fprintf(fid3,fmt1e,' '); % Blank line
fprintf(fid3,fmt1a,'Ave     ');
fprintf(fid3,fmt1d,astats(1,:));



% The unadjusted reference-period monthly standard devs for stations 
% and the global average standard deviation over stations

file3=['ustdev' strsuf '.txt'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt1a='%s ';
fmt1e='%s\n';
fmt1b='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f';
fmt1c='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n';
fmt1d=[fmt1b fmt1c];
	
for n=1:ns2;
	b3=B3(n,:);
	x=ustdev(n,:);
	fprintf(fid3,fmt1a,b3);
	fprintf(fid3,fmt1d,x);
end
fprintf(fid3,fmt1e,' '); % Blank line
fprintf(fid3,fmt1a,'Ave     ');
fprintf(fid3,fmt1d,ustats(2,:));

fclose(fid3);


% The adjusted reference-period monthly standard devs for stations 
% and the global average standard deviation over stations

file3=['astdev' strsuf '.txt'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt1a='%s ';
fmt1e='%s\n';
fmt1b='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f';
fmt1c='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n';
fmt1d=[fmt1b fmt1c];
	
for n=1:ns2;
	b3=B3(n,:);
	x=astdev(n,:);
	fprintf(fid3,fmt1a,b3);
	fprintf(fid3,fmt1d,x);
end
fprintf(fid3,fmt1e,' '); % Blank line
fprintf(fid3,fmt1a,'Ave     ');
fprintf(fid3,fmt1d,astats(2,:));

fclose(fid3);



% Sample size (number of years of valid data in reference period)
% for computing means and standard deviations

file3=['usize' strsuf '.txt'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt1a='%s ';
fmt1e='%s\n';
fmt1b='%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f';
fmt1c='%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n';
fmt1d=[fmt1b fmt1c];
	
for n=1:ns2;
	b3=B3(n,:);
	x=usize(n,:);
	fprintf(fid3,fmt1a,b3);
	fprintf(fid3,fmt1d,x);
end

fclose(fid3);


% The ratios of monthly means used to adjust station means 
file3=['ratmn' strsuf '.txt'];
pf3=[file3];
fid3=fopen(pf3,'w');

fmt1b='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f';
fmt1c='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n';
fmt1d=[fmt1b fmt1c];
	
for n=1:ns2;
	b3=B3(n,:);
	x=ratmn(n,:);
	fprintf(fid3,fmt1a,b3);
	fprintf(fid3,fmt1d,x);
end

fclose(fid3);


% The ratios of monthly std devs used to adjust station std devs 
file3=['ratsd' strsuf '.txt'];
pf3=[file3];
fid3=fopen(pf3,'w');

fmt1b='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f';
fmt1c='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n';
fmt1d=[fmt1b fmt1c];
	
for n=1:ns2;
	b3=B3(n,:);
	x=ratsd(n,:);
	fprintf(fid3,fmt1a,b3);
	fprintf(fid3,fmt1d,x);
end

fclose(fid3);

% Time series matrix of number of stations in the regional average
file3=['nstns' strsuf '.dat'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt2a='%4.0f';
fmt2e='%s\n';
fmt2b='%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f';
fmt2c='%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n';
fmt2d=[fmt2b fmt2c];
	
for n=1:nyrs3;
	year=yr3(n);
	x=NU(n,:);
	fprintf(fid3,fmt2a,year);
	fprintf(fid3,fmt2d,x);
end

fclose(fid3);


% Time series matrix of regional stdzd anomalies, unadjusted
file3=['uregsn' strsuf '.dat'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt2a='%4.0f';
fmt2e='%s\n';
fmt2b='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f';
fmt2c='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n';
fmt2d=[fmt2b fmt2c];
	
for n=1:nyrs3;
	year=yr3(n);
	x=ZU(n,:);
	fprintf(fid3,fmt2a,year);
	fprintf(fid3,fmt2d,x);
end

fclose(fid3);


% Time series matrix of regional stdzd anomalies, adjusted
file3=['aregsn' strsuf '.dat'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt2a='%4.0f';
fmt2e='%s\n';
fmt2b='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f';
fmt2c='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n';
fmt2d=[fmt2b fmt2c];
	
for n=1:nyrs3;
	year=yr3(n);
	x=ZA(n,:);
	fprintf(fid3,fmt2a,year);
	fprintf(fid3,fmt2d,x);
end

fclose(fid3);


% Time series matrix of regional climatic series, converted from
% stdzd anomalies back to original units: unadjusted
file3=['ureg' strsuf '.dat'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt2a='%4.0f';
fmt2e='%s\n';
fmt2b='%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f';
fmt2c='%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n';
fmt2d=[fmt2b fmt2c];
	
for n=1:nyrs3;
	year=yr3(n);
	x=YU(n,:);
	fprintf(fid3,fmt2a,year);
	fprintf(fid3,fmt2d,x);
end

fclose(fid3);


% Time series matrix of regional stdzd anomalies, adjusted
file3=['areg' strsuf '.dat'];
pf3=[file3];
fid3=fopen(pf3,'w');

blnk=' ';
fmt2a='%4.0f';
fmt2e='%s\n';
fmt2b='%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f';
fmt2c='%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n';
fmt2d=[fmt2b fmt2c];
	
for n=1:nyrs3;
	year=yr3(n);
	x=YA(n,:);
	fprintf(fid3,fmt2a,year);
	fprintf(fid3,fmt2d,x);
end

fclose(fid3);


% .mat files of the regional climatic series (same data as ureg*.dat and
% areg*.dat
flout=['ureg' strsuf];
WU=[yr3 YU];
WA=[yr3 YA];
eval(['save ' flout ' WU']);
flout=['areg' strsuf];
eval(['save ' flout ' WA']);



fclose all;


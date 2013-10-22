function regcli3(flsin,miscin,Sinf,A,pntnum)
%
% Standardized anomaly method for regional climate series -- point-by-point
% version.  Adapted from regcli2.m  Regcli3.m interpolates station PPT or TMP 
% monthly data to gridpoints, or to any network of lon-lat points.  Main
% differences from regcli2.m:
%
%	-regcli3 has no user prompted info; info is passed as arguments
%	-regcli3 does not implement adjustment of reference period means or standard
%		devs. In terminology of regcli2.m, a 'pass-1' analysis is all that is
%		done
%	-regcli3 is intended to be called for individual gridpoints (or, say, tree
%		tree-ring sites) for which interpolated series are desired.  For a 
%		western North America application, calling function would call regcli3.m
%		repeatedly to build up interpolated series
%
% Application.  Wanted tree-site-centered TMP series at 700+ western NA tree
% sites, from network of 1002 GHCN v2 TMP stations. Direct use of regcli2.m 
% for this problem would be difficult because of huge array sizes. For example,
% 1002 x 700 just to store the station-to-tree distances, and 
% 100 yr x 700 stations x 12 months to store all stations' monthly data.
% Decided instead to write regcli3.m to treat problem tree-site by tree-site
%
% Weighting implemented so far is  inverse distance
% 
% Meko 11-16-97
%
%
%**************************** IN ARGS *****************************
%
% flsin (? x 5?)s filenames (w/o .mat) of .mat files with monthly clim
% miscin {} miscellaneous input
%		{1} CG [1 x 2]r long-lat coords of gridpoint.
%		{2} C2 [ns1 x 2]r long-lat file in mapping units for the ns1 stations
%			from which gridpoint data is to be interpolated
%		{3} A (1 x ?)r  distance (km) from target point to the ns1 stations
%		{4} yrinf (2 x 2)i  First and last year of :
%			row 1: reference period (start, end years) for computing
%					monthly means and standard deviations to be used in computing
%					standardized anomalies
%			row 2: storage matrix to hold monthly data for all stations; this
%					period  must cover at least the period with any valid data
%					at any station
%		{5} kopts (1 x 4)i options
%			kopts(1): mode for running program
%					==1: diagnostic (lots of feedback)
%					==2: operational (just the facts, Ma'am)
%			kopts(2):==1 --> subst. zero for negative monthly regional values
%				    ==2 --> accept negative monthly regional values
%					(note: typically will set ==1 for pcp, ==2 for tmp)
%			kopts(3) ==1 nearest nmax stns, disregarding search radius is to
%							be used in weighting gridpoints.  Will always have
%							nmax stations in weighting, unless fewer than nmax
%							are available in that month/year for the whole ns1
%							station set
%				    	==2 nearest up-to-nmax stns,  and all must be in
%							search radius;  might need few trial and error
%							runs with different dsrch values to make sure
%							all gridpoints can be estimated
%			kopts(4) ==1 weighting method inverse distance
%							(only option so far)
%
%		{6} (1 x 3)r distance settings	
%			dcrit -- critical threshold distance to avoid over-weighting site in
%				distance weighting (km)
%			dsrch -- search distance for weighting (km)
%			nmax --  max number of stations to weight to a gridpoint
%		{7} paths{} string variables with paths
%			(1) to input files of monthly climate data
%			(2) to output files  of interpolated monthly gridpoint data
%
%	Sinf {} string matrices of identifying information
%		{1} (1 x 70)s for the target point
%		{2} (? x 70)s for the climate stations represented by the .mat files
%
%	A (1 x ?)r  distance (km) for target point to each station (ns1 stations)
%	pntnum(1 x 1)i  point number -- numeric identifier used to build prefixes
%			of output files
%
%**************************** OUTPUT FILES ***********************************
%
% A .mat file holding the target point monthly climate series in matrix Z
% The name of this .mat file is built from the input variable pntnum, as
% ['P' int2str(pntnum)] for pcp, and ['T' int2str(pntnum) for tmp
%
% An ascii information file for the point, with same prefix as the .mat file,
% but with .txt extension.  The first line of this file contains information
% on the target point itself: name of .mat data file, tree-ring site name (if
% a tree-ring site), lon and lat, elevation in m, first and last year with
% any monthly data.  Lines 2-6 contain information on the climate stations used
% to form the target-point reference-period means and standard deviations.
% These are usually the 5 nearest stations, or sometimes fewer than 5 depending
% on kopts(3). Distance alone is not the only criterion.  The station record
% must also have the minimum required percentage of valid data for each month
% in the reference period. Information on lines 2-6 is seq number (nearest
% as #1); Station .mat file name, 
% station name, state, country, lon/lat, elevation, time covrage
% distance to target point, and weight for the target-point mean and std dev Example:
% 
%	p1 		Tsegi Point PIED  	-110.32	34.56  1243m 1904 1995
%	1-D67		Holbrook AZ,USA    	-109.23 35.64   934m 1912 1996 32km  0.123
%
%  
% INDIVIDUAL STATION INFO
%
% stnlst?.txt --   seq number; 4-letter code; long, lat;
%		(i) if in region (i), () if not used; "*" if a master series
% mnstn?.txt --  stn code; 12 reference-period monthly means
% sdstn?.txt --  stn code; 12 ref-pd standard devs
% yrsref?.txt -- stn code; 12 sample sizes -- number of years ref-pd means
%			and standard deviations based on
% ratmn?.txt -- stn code; 12 ratios -- ratio of full ref pd mean to 
%		subperiod mean in the master series. These ratios may or may not have
%		been used to adjust the station means in computing stdzd anomalies, 
%		depending on kopts.
% ratsd?.txt -- ttn code; 12 values -- like ratmn?.txt, but for standard devs
%
%
% GRIDPOINT INFO
%
% Files Pnnnn.???, where nnnn is the  
% grdlst?.txt --  gridpoint #; long, lat; (j); weight on for regioinal series
% mngrd?.txt --  gridpoint #, 12 reference-period monthly means
%		These means by distance weighting station means to gridpoints as
%		listed in wgts?.txt
% sdgrd?.txt --  gridpoint #; 12 ref-pd standard devs, weighted as in mngrd?.txt
% wgts?.txt -- How station means and std devs weighted to gridpoint. 
%		For each gridpoint, line 1 lists the point # and long/lat
%		Succeeding lines list  stn sequence no., station code, and weight
% grdm?j.mat -- weighted climate monthly time series for gridpoint j
%		Data in Z, with year as column 1
% grdm?.mat -- multi-gridpoint file of weighted climate monthly time series for 
%			gridpoints.  Data in X.  Year as col 1.  Cols 2-13 apply to point1,
%			cols 14-25 to point 2, etc.
%
%
% REGION INFO
%
% mnreg?.txt --  region #, 12 reference-period monthly means
%		These means by latitude-weighting gridpointon means as
%		listed in grdlst?.txt
% sdreg?.txt --  region #; 12 ref-pd standard devs, weighted as grdlst?.txt
% nstns?i.dat -- each row a year and 12 sample sizes (for months of year)
%		telling how many stations in the arbitrarily specified region i in 
%		each month/year.  Station-region makeup specified in input IR1.
% nstns?.dat -- combined information on number of stations, combined over all
%		regions
% mstsa?.mat or regsa?.mat -- multi-region monthly standardized anomaly matrix.
%		Col 1 is year. Next 12 cols are Jan-Dec for region 1; next 12 for region
%		2, ...  File is "mst" if pass 1 according to kpass. File is "reg" if 
%		pass 2.  Data in matrix X.
% mst?.mat  or reg?.mat -- like mstsa?.mat/regsa?.mat, except transformed back
%		to climatic units using the weighted means and standard deviations in
%		mnreg?.txt and sdreg?.txt. Data in matrix X.
% mst?i.mat  -- individual 13-col files of regional climate series extracted
%		from mst?.mat.  Data in Z.
%
%
%************************* USER WRITTEN FUNCTIONS CALLED ***************
% gcdist.m great circle distance
% grd2reg.m gridpoint data latitude-weighted to regional
% regsub1.m  compute number of stations matrices
% wgtdist1.m inverse distance weighting
%
%******************************** NOTES *************************************
%
%
% Reference: Jones & Hulme (1996).  
% Method.   Monthly station time series are first converted to standardized
% anomaly series using the mean and standard deviation for that station/month
% computed for a specified common reference period.  NaN data are ignored in
% computing the means and standard deviations.
%
% The standardized anomaly series are eighted over stations around the taget 
% point  (subject to constraint on search radius) to produce a single
% interpolated average standardized anomaly monthly time series.  The weighting 
% function may vary year-to-year depending on stations available
%
% The reference-period station means and standard deviations are similarly
% weighted to the target point.  The nmax nearest stations are used in the 
% weighting if kopt(3)==1, or the nearest up-to-nmax stations within dsrch of
% the target point if kopt(3)==2
% 
% The target point standardized anomalies are converted to original units
% (e.g., inches for ppt) using the interpolated target point monthly
% reference period means and standard deviations
%
% The "target point" might be a desired long-lat gridpoint if the objective is
% to grid station climate data to a grid, or might be a tree-ring site if the
% goal is a "tree-centered" climate series
%
%******************** DEFINTIONS OF COMPUTED VARIABLES *******************
%
% A distance mtx, gridpoints to stations (km)
% B1 (ns1 x 8)s .mat file prefixes
% D huge monthly mtx, original data
% E huge ..., unadjusted standardized anomalies
% mstpref... prefix of master series .mat files
% ncols -- number of cols in huge monthly multi-station matrix
% npts1 (rv)  -- number of gridpoints in each region
%- nmast number of stations in master network
% nref -- # yr in reference period
%- nreg -- number of regions
% nrows -- number of rows in huge monthly multi-station matrix
%- ns1 - total number of .mat monthly input series
% ns2 (rv) number of stations in each region
%- path3 path to raw monthly climate .mat files
% path5 path to output files made by regcli2.m
% P (nrows x (mG*12)) huge mtx of gridpoint monthly series, orig units
% PM (mG x 12) long-term means at gridpoints
% PS (mG x 12) long-term std devs at gridpoints
% path3 -- path to input monthly .mat files
% Q like P, except stdzd anomaly version
% TMN,TSD -- stores AMN,ASD or UMN,USD depending on kopts(5)
% UMN unadjusted long-term means, all stations
% USD adjusted long-term means all stations
% wt2 -- weights on gridpoints to regions
% W5 -- weights on stations to gridpoints for long-term means
% yr1 -- year vector for master reference period
% yr2 -- year vector for huge multistation monthly matrix


%
% Allocate
blnks=blanks(8);
B1=blnks(ones(200,1),:); % Initialize to allow for 200 filenames
sblank=sprintf('%\n'); % used for blank line to screen



% Prompt for whether pass1 or pass 2.
filecode=input('One-letter filecode for record keeping: ','s');


%********************************************************************
% Read list of station .mat files; count # of files; store file names
%First line of the list is the path to input files
%Second line is the path to output files
[file1,path1]=uigetfile('fls*.txt','List of input .mat files');
pf1=[path1 file1];
fid1=fopen(pf1,'r');
path3=strtok(fgetl(fid1));
path5=strtok(fgetl(fid1));
k1=1;
ns1=0; % counter for number of .mat files in list
while k1;
	c1=fgetl(fid1);
	if feof(fid1) | length(c1)<4;
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
clear k1 c1 c2 nc2 
fclose(fid1);

% B1 now contains file name prefixes
% ns1 is the number of files
% path3 is the path to the .mat input station data files


%*********** QC the input args

[file2,path2]=uigetfile('misc*.mat','Miscellaneous program control');
pf2=[path2 file2];
eval(['load ' pf2]);

% Check that all variables wanted from the misc file exist
Ltemp1=[exist('CG') exist('C2')];
Ltemp2=[exist('pdref') exist('pdstore')];
Ltemp3=[exist('kmode') exist('kopts')];
Ltemp=[Ltemp1 Ltemp2 Ltemp3];
Ltemp=all(Ltemp==1);
if ~Ltemp
	clc
	disp('The following variables should be in the misc.mat file:');
	disp('CG,C2,pdref,pdstore,kmode,kopts');
	error('All the above were not in the file');
end
clear Ltemp1 Ltemp2 Ltemp3 Ltemp file2

% Check input arguments


[mG,nG]=size(CG); % x,y coordinates of target point
if nG~=2 | mG~=1 ;
	error('CG must be 1 x 2')
end
if  any(any (CG(:,1)>0));
	error('CG long should be negative')
end
if  any(any(CG(:,2)>90));
	error('CG(:,2) should be all 90 deg or less (latitude)');
end
clear nG

[mC2,nC2]=size(C2);
if mC2~=ns1;
	error('Row size of C2 not equal to number of .mat filenames');
%end
if nC2~=2;
	error('C2 should be 2-column');
end
if  any(any (C2(:,1)>0));
	error('C2 long should be negative')
end
if  any(any(C2(:,2)>90));
	error('C2(:,2) should be all 90 deg or less (latitude)');
end
clear mC2,nC2


npts = 1;  % Number of 'target' points


% pdref holds start, end years for reference period
if pdref(1)>pdref(2); 
	error('Reference period ends before it starts');
end
yr1=(pdref(1):pdref(2))';
nref=length(yr1); % # yrs in ref pd

% pdstore holds start, end years of a matrix to store the monthly
% data for all stations. 
if pdstore(1)>pdstore(2)
	error('end year of pdstore precedes start year');
end

% Mode is either rapid or diagnostic
if kmode<1 | kmode>2,
	error('Invalid setting for kmode')
end

% kopts holds program options
if ~isnan(kopts(1));
	error('kopts(1) should be set to NaN');
if kopts(2)<1 | kopts(2)>2;
	error('Illegal setting of kopts(2)')
end
if kopts(3)<1 | kopts(3)>2;
	error('Illegal setting of kopts(3)')
end
if kopts(4)<1 | kopts(4)>1;
	error('Illegal setting of kopts(4)')
end

regpref='mst'; % prefix to output "regional" .mat files


disp(' ');
disp('Press any key to continue');
pause

clc



% Optional feedback on station composition in regions
if kmode==1; % diagnostics mode
	clc
	disp('Listing of all monthly climate stations')
	
	disp(B1);
	
	disp('Press any key to continue')
   pause
   clc
	
	disp('Data:');
	disp([int2str(ns1) ' stations']);
	disp([int2str(nreg) ' regions']);
	disp([int2str(mG) ' gridpoints'])
	
	disp('Number of stations:')
	disp(ns2)
	
	disp('Reference period')
	disp([pdref(1) pdref(2)]);
	
	disp('Period for huge matrix');
	disp([pdstore(1) pdstore(2)]);
	
	disp('Settings for call to wgtdista.m')
	disp([' dcrit = ' num2str(dcrit) ' km'])
	disp([' dsrch= ' num2str(dsrch) ' km'])
	disp([' nmax= ' int2str(nmax) ' stations'])
	disp([' Option for ignoring search radius: ' int2str(kopts(3))])
	disp([' Option for weighting method: ' int2str(kopts(4))])
	 
	disp('Press any key to continue')
   pause
   clc
end

% Allocate
a=NaN;
UMN=a(ones(ns1,1),ones(12,1));
USD=a(ones(ns1,1),ones(12,1));
USIZE=a(ones(ns1,1),ones(12,1));

%**************** READ IN MONTHLY STATION CLIMATE SERIES
%
disp('READING MONTHLY DATA FILES AND STORING DATA');
disp(' ');
% Store the input monthly ppt series for all specified climate stations for
% the region in a single matrix

% Initialize matrices 

yr2=(pdstore(1):pdstore(2))'; % years vector for storage matrix
nrows = length(yr2); % number of rows (years) in storage matrix
ncols = ns1*12; % % number of cols in storage matrix
D=a(ones(nrows,1),ones(ncols,1));  % monthly raw data
E=a(ones(nrows,1),ones(ncols,1));  % monthly data, unadjusted anoms
T=a(ones(ns1,1),ones(2,1));  % start, end years of each station's record
I=a(ones(ns1,1),ones(2,1));  % corresp row pointer to D


% What cols in D and E will the station monthly data go into?
j1=(1:12:(ncols-11)); % start col for stn1, 2,etc
j2=j1+11; % end col
J1=[j1' j2'];
% First row of J1 will give first and last storage row for stn 1, second row
%		stn 2 etc

% What cols in D and E will hold all Jan data, Feb data, etc?
j1=(1:12);
j1=j1(ones(ns1,1),:);
j2=(1:ns1)';
j2=12*(j2-1);
j2=j2(:,ones(12,1));
J2=j1+j2;
% First col of J2 will point to cols of D for Jan, secnd for Feb, etc


% Loop over stations: get first and last year,
% 
for n=1:ns1;
	clear X  Y  Z;
	fln=strtok(B1(n,:));
	eval(['load ' path3 fln]);
	if (exist('Z')~=1  & exist('X')~=1  & exist('Y')~=1)
		error(['No variable named X, Y or Z in: ' fln]);
	elseif exist('X')==1 & exist('Y')~=1  & exist('Z')~=1; % X is the variable
		Z=X;  
		clear X;
	elseif exist('Y')==1 & exist('X')~=1  & exist('Z')~=1; % Y is the variable
		Z=Y;  
		clear Y;
	elseif exist('Z')==1 & exist('X')~=1  & exist('Y')~=1; % Z is the variable
	else
		error('data file holds more than one of X, Y or Z');
	end
	

	% Make year vector; store first year and last year; calc and store
	% row reference info
	yr = Z(:,1);
	nyrs=length(yr);
	if yr(1)<pdstore(1) | yr(nyrs)>pdstore(2)
		disp(fln)
		error('Above file has data outside period specified as pdstore');
	end
	T(n,:)=[yr(1) yr(nyrs)]; % first and last years of station data
	I(n,:)=[yr(1)-pdstore(1)+1  yr(nyrs)-pdstore(1)+1]; % row index into D and E


	% Store the monthly data
	D(I(n,1):I(n,2),J1(n,1):J1(n,2)) = Z(:,2:13);
end


%	*******************  MEAT OF THE ORDER  *************************************
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

umean=a(:,ones(12,1));  % to hold reference period station means for a station
ustdev=a(:,ones(12,1)); % ... standard deviations
usize=a(:,ones(12,1)); % ... sample size

igo=pdref(1)-pdstore(1)+1;
isp=pdref(2)-pdstore(1)+1;
irows=(igo:isp);

for nn=1:ns1; % loop over all stations in the region
	jthis=nn; % index to station number
	jcols=J1(jthis,1):J1(jthis,2);
	X=	D(irows,jcols); % reference period monthly data
	L1=~isnan(X);
	usize(nn,:)=sum(L1);  % store sample size (yrs in ref pd with valid data)
	[mn1,msize]=meannan(X); % unadjusted ref pd means; second arg not used
	[std1,msize]=stdnan(X); % unadjusted std devs
	umean(nn,:)=mn1;
	ustdev(nn,:)=std1;

	% Store reference period means and standard deviations for this station
	UMN(jthis,:)=umean(nn,:);
	USD(jthis,:)=ustdev(nn,:);
	USIZE(jthis,:)=usize(nn,:);
end


%********************************************************************
%
% Compute and store means and standard deviations averaged over
% all stations These arithmetic means are not used elsewhere in the program, but 
% might be interesting to compare with the target-point monthly means
% and standard devs by the distance weighting algorithm
ustats=a(ones(2,1),ones(12,1)); % storage for means in row 1, sdevs in 2
ustats(1,:)=mean(UMN);
ustats(2,:)=mean(USD);


%*******************************************************************
%
% Compute standardized departure time series for stations.
% Fill matrices E and F with standardized departure series. 

for nn=1:ns1; % Loop over stations
	ithis=nn;
	irows=I(ithis,1):I(ithis,2); % points to rows in D holding this station's data
	jcols=J1(ithis,1):J1(ithis,2); % ...cols in D ....
	X=D(irows,jcols); % Monthly data, entire data record for this station
	[mX,nX]=size(X);

	mns = UMN(nn,:); % monthly ref-period means for this station
	sds = USD(nn,:); % ... standard devs ...
	MNS=mns(ones(mX,1),:); % rv to matrix
	SDS=sds(ones(mX,1),:); % ditto
	SA = (X - MNS) ./ SDS; % standardized anomalies for this station
	E(irows,jcols)= SA;  % store with standzd anomalies of other stations
end



%***************  COMPUTE TARGET POINT WEIGHTED STDZD ANOMALIES

disp('WEIGHTING STDZD ANOMALIES TO TARGET POINT');

% Allocate
P=a(ones(nrows,1),ones(12,1)); % to hold target point data in clim units
Q=a(ones(nrows,1),ones(12,1)); % to hold target point data in stdzd units

% Build pointers for each gridpoint's series in P and Q

% What cols in P and Q will the target point's monthly data go into?
JP1=[1 12]; % first and last cols

% What cols in P and Q will hold all Jan data, Feb data, etc?
JP2= [1:12]; % says Jan data in col 1, Feb in col 2, etc
% First col of JP2 will point to cols of P with Jan data , 
% second col to cols of P with Feb data, etc

% Recall that
%
%	A is a rv of distances (km) from the target point to the ns1 stations

% Loop over the 12 months of the year
for kmon=1:12;

	disp(['   Month ' int2str(kmon)]);
	jp2=JP2(:,kmon); %cv of cols indices to P,Q for this month

	% Pull the subset of station stdzd anomalies, all yr in E 
	F1=E(:,J2(:,kmon));
	

	% Pull subset of rows of F1 for which not all values are NaN
	L1=~isnan(F1);
	L1a=(any(L1'))';
	F2=F1(L1a,:);

	% Make logical pointer to F2 with data
	L1=~isnan(F2);

	% Compute matrix of weights for this month of year
	kk=[kopts(3) kopts(4)];
	[WW,NW,JW]=wgtdist1(A,L1,dcrit,dsrch,nmax,kk);

	if kmode==1;
		disp(['GRIDPOINT # ' int2str(kgrd)]);
	end

	jpthis=jp2(kgrd); % col in P, Q for this month		

	colsWW=JW(kgrd,1):JW(kgrd,2);  % columns of W with this points wts
	W1=WW(:,colsWW);
		
	% Weight the stations
	G1=F2 .* W1;
      
	% Sum over cols to get cv of weighted stdzd anoms this point/month
	g1 = (sumnan(G1'))';

	% Store weighted anomalies
	Q(L1a,jpthis)=g1;
end



%************** GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS
%
disp('COMPUTING GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS');

% Allocate
PM=a(:,ones(12,1)); % to hold means
PS=a(:,1),ones(12,1)); % to hold std devs
W5=a(:,ones(ns1,1));

TMN=UMN;  % means
TSD=USD; % standard devs


% Mean
L1=~isnan(TMN');
[WWM,NWM,JWM]=wgtdist1(A,L1,dcrit,dsrch,nmax,kk);
% Note that all rows of WWM are identical. Only need first
for n=1:mG;
	jcols=JWM(n,1):JWM(n,2);
	W5(n,:)=WWM(1,jcols);
end
L1=isnan(W5);
sum5=sum(sum(L1));
% Matrix mult with NaNs gives Nans, so change NaNs to zero
if sum5>0;
	W5(L1)=zeros(sum5,1);
end

% Make temporary matrices TTMN and TTSD by replacing NaNs with zeros
% in TMN, TSD before matrix mult.  Otherwise, NaNs in matrix mult, which
% means all products NaN
TTMN=TMN;
TTSD=TSD;
L1=isnan(TMN);
sum5=sum(sum(L1));
if sum5>0;
	TTMN(L1)=zeros(sum5,1);
end
L1=isnan(TSD);
sum5=sum(sum(L1));
if sum5>0;
	TTSD(L1)=zeros(sum5,1);
end

PM=W5*TTMN;
PS=W5*TTSD;
clear WWM NWM JWM L1 sum5 TTMN TTSD


%********* CONVERT GRIDPOINT STANDARDIZED ANOMALIES TO ORIGINAL UNITS

[mQ,nQ]=size(Q);
PM1=repmat(PM,mQ,1);
PS1=repmat(PS,mQ,1);
P = (Q .* PS1) + PM1;

% Optionally substitute zero for negative monthly regional values
if kopts(2)==1; % substitute -- most appropriate for ppt data
	L3=P<=0;
	numneg=sum(sum(L3));
	if numneg>0;
		P(L3)=zeros(numneg,1);
	end
elseif kopts(2)==2; % do not substitute -- most approp. for tmp data
	% no action
else
	error('Invalid setting for kopts(2)');
end


 
%*********** COUNT NUMBER OF STATION WITH DATA IN EACH MONTH/YEAR  ***********
disp('CALLING regsub1 TO GET STATION COUNT MATRIX  NS')
[NS,yr3]= regsub2(ns1,E,J2,yr2);


%******************** *****************************
disp('BUILDING OUTPUT FILES')
 
%****** OUTPUT THE RESULTS **************************
%
% 
% RECALL THAT
%
% B1 station filename list
% C2 long lats of stations
% CG long lats of gridpoints
% NS sample size for target point series for each month/year
% P tsm of monthly target point data, original units
% Q tsm ... as standardized anomalies
% PM long term monthly means of target point data
% PS long term standard devs
% TMN long-term (ref period) monthly means of station data
% TSD long-term   std devs
% W5 matrix of weights used to compute long-term target point
%	  means from long term station means
% yr2 year vector for huge matrix
% yr3 year vector for NS

% filecode -- 1-char letter added to file prefixes for record keeping



% Reduce row size of P, Q, by lopping off years for which
% all elements are NaN
L1=isnan(Q);
L2=(~all(L1'))';
yr4=yr2(L2);
nyrs4=length(yr4);
P1=P(L2,:);
Q1=Q(L2,:);

L1=isnan(RA);
L2=(~all(L1'))';
yr7=yr2(L2);
nyrs7=length(yr7);
RA1=RA(L2,:);
R1=R(L2,:);


% Station list, with sequence number, long, lat, region station in
% for adjustment of mean and standard dev, and number indicating
% if station a master -- and if so, for what region
file3=['stnlst' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
blnk='  ';
for n=1:ns1;
	% get region number
	[i,j]=find(IR1==n); % j is region
	str1=sprintf('%3.0f  ',n); % station number
	str2=sprintf('%s',B1(n,:)); % station code (file prefix)
	str3=sprintf('%6.2f %5.2f ',C2(n,[1 2])); % long, lat
	str4=sprintf('%2.0f',j); % region number
	% master or not
	[i,j]=find(IM==n);
	if isempty(j); % not a master
		str5=sprintf('%s',blnk);
	else
		str5=sprintf('%s','*');
		%str5=sprintf('%2.0f',j); % used to want region in parens
	end
	strall=[str1 str2 str3 '(' str4 ')' str5];
	fprintf(fid3,'%s\n',strall);
end
fclose(fid3);


% Gridpoint list:
% gridpoint #; region gridpoint in; weight of gridpoint in region
file3=['grdlst' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
blnk='  ';
for n=1:mG;
	% get region number
	[i,j]=find(IR2==n); % j is region
	str1=sprintf('%3.0f  ',n); % gridpoint number
	str2=sprintf('%6.2f %5.2f ',CG(n,[1 2])); % long, lat
	str3=sprintf('%2.0f',j); % region number
	% get weight
	w6=wt2(n);
	str4=sprintf('%6.4f',w6); % weight
	strall=[str1 str2 '(' str3 ') ' str4];
	fprintf(fid3,'%s\n',strall);
end
fclose(fid3);

% Listing of stations and weights used to get long term
% gridpoint means and std devs
file3=['wgts' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
blnk='  ';
for n=1:mG;
	str1=sprintf('%3.0f  ',n); % gridpoint number
	str2=sprintf('%6.2f %5.2f ',CG(n,[1 2])); % long, lat
	stra=[str1 str2];
	fprintf(fid3,'%s\n',stra);
	% get the weights vector and find number of weights
	w5=W5(n,:); % rv, weights 
	% find the non-zero elements
	j=find(w5~=0);
	w6=w5(j);
	nwgts=length(w6); % # weights
	for k = 1:nwgts;
		str3=sprintf('  %3.0f',j(k));
		str4=sprintf('  %s',B1(j(k),1:5));
		str5=sprintf('(%7.4f)',w6(k));
		strb=[str3 str4 str5];
		fprintf(fid3,'%s\n',strb);
	end
end
fclose(fid3);


% .mat file of multi-region monthly stdzd anomalies
flout=[regpref 'sa' filecode];
X=[yr7 RA1];
eval(['save ' path5 flout ' X']);

% .mat file of multi-region  monthly series in original units
flout=[regpref filecode];
X=[yr7 R1];
eval(['save ' path5 flout ' X']);


% .mat file of multi-gridpoint  monthly series in original units
flout=['grd' filecode];
X=[yr4 P1];
eval(['save ' path5 flout ' X']);


% individual .mat files of monthly regional climate series
for n=1:nreg;
	% build file name
	flout=[regpref filecode int2str(n)];
	% get cols of data for region
	jcols=JP3(n,1):JP3(n,2);
	Z=R(:,jcols);
	% truncate any rows all NaN
	L1= isnan(Z);
	L2=   (~all(L1'))';
	yr=yr2(L2);
	Z=Z(L2,:);
	% Put year as col 1
	Z=[yr Z];
	% output
	eval(['save ' path5 flout ' Z']);
end


% Build and save individual gridpoint .mat and .dat files
for n = 1:mG; % loop over gridpoints
	% build file name
	flout=['grd' filecode int2str(n)];
	% get cols of data for gridpoint
	jcols=JP1(n,1):JP1(n,2);
	Z=P(:,jcols);
	% truncate any rows all NaN
	L1= isnan(Z);
	L2=   (~all(L1'))';
	yr=yr2(L2);
	Z=Z(L2,:);
	% Put year as col 1
	Z=[yr Z];
	% output
	eval(['save ' path5 flout ' Z']);
end

% Ref period station means.  These might be adjusted or
% unadjusted, depending on kopts(1) and kopts(5) settings
% Note that on pass 1, all except the master series would be
% NaN.  Want to cull out the non-NaN rows
X=TMN;
L1=isnan(X);
L2=(~all(L1'))';
X=X(L2,:);
BB1=B1(L2,:);
mX=size(X,1);
file3=['mnstn' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
if kopts(5)==1;
	fprintf(fid3,'%s\n\n','These means are adjusted');
else
	fprintf(fid3,'%s\n\n','These means are unadjusted');
end
blnk=' ';
fmt1a='%s';
fmt1e='%s\n';
fmt1b='%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f';
fmt1c='%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n';
fmt1d=[fmt1b fmt1c];
for n=1:mX;
	x=X(n,:);
	fprintf(fid3,fmt1a,BB1(n,1:5));
	fprintf(fid3,fmt1d,x);
end
fclose (fid3);

% Ref period gridpoint monthly means .  These might be adjusted or
% unadjusted, depending on kopts(1) and kopts(5) settings
file3=['mngrd' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
if kopts(5)==1;
	fprintf(fid3,'%s\n\n','These means are adjusted');
else
	fprintf(fid3,'%s\n\n','These means are unadjusted');
end
blnk=' ';
fmt1a='%s';
fmt1e='%s\n';
fmt1b='%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f';
fmt1c='%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n';
fmt1d=[fmt1b fmt1c];
for n=1:mG;
	x=PM(n,:);
	fprintf(fid3,'#%3.0f',n);
	fprintf(fid3,fmt1d,x);
end
fclose(fid3);

% Ref period region monthly means .  These might be adjusted or
% unadjusted, depending on kopts(1) and kopts(5) settings
file3=['mnreg' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
if kopts(5)==1;
	fprintf(fid3,'%s\n\n','These means are adjusted');
else
	fprintf(fid3,'%s\n\n','These means are unadjusted');
end
blnk=' ';
fmt1a='%s';
fmt1e='%s\n';
fmt1b='%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f';
fmt1c='%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n';
fmt1d=[fmt1b fmt1c];
for n=1:nreg;
	x=RM(n,:);
	fprintf(fid3,'#%3.0f',n);
	fprintf(fid3,fmt1d,x);
end
fclose(fid3);

% Ref period station std devs.  These might be adjusted or
% unadjusted, depending on kopts(1) and kopts(5) settings
file3=['sdstn' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
if kopts(5)==1 kopts(1)==2;;
	fprintf(fid3,'%s\n\n','These standard devs are adjusted');
else
	fprintf(fid3,'%s\n\n','These standard devs are unadjusted');
end
blnk=' ';
fmt1a='%s';
fmt1e='%s\n';
fmt1b='%6.3f%6.3f%6.3f%6.3f%6.3f%6.3f';
fmt1c='%6.3f%6.3f%6.3f%6.3f%6.3f%6.3f\n';
fmt1d=[fmt1b fmt1c];
% cull non-NaN rows
X=TSD;
L1=isnan(X);
L2=(~all(L1'))';
X=X(L2,:);
mX=size(X,1);
BB1=B1(L2,:);
for n=1:mX;
	x=X(n,:);
	fprintf(fid3,fmt1a,BB1(n,1:5));
	fprintf(fid3,fmt1d,x);
end
fclose(fid3);

% Ref period gridpoint std devs.  These might be adjusted or
% unadjusted, depending on kopts(1) and kopts(5) settings
file3=['sdgrd' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
if kopts(5)==1 kopts(1)==2;;
	fprintf(fid3,'%s\n\n','These standard devs are adjusted');
else
	fprintf(fid3,'%s\n\n','These standard devs are unadjusted');
end
blnk=' ';
fmt1a='%s';
fmt1e='%s\n';
fmt1b='%6.3f%6.3f%6.3f%6.3f%6.3f%6.3f';
fmt1c='%6.3f%6.3f%6.3f%6.3f%6.3f%6.3f\n';
fmt1d=[fmt1b fmt1c];
for n=1:mG;
	x=PS(n,:);
	fprintf(fid3,'#%3.0f',n);
	fprintf(fid3,fmt1d,x);
end
fclose(fid3);


% Ref period region std devs .  These might be adjusted or
% unadjusted, depending on kopts(1) and kopts(5) settings
file3=['sdreg' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
if kopts(5)==1 kopts(1)==2;;
	fprintf(fid3,'%s\n\n','These standard devs are adjusted');
else
	fprintf(fid3,'%s\n\n','These standard devs are unadjusted');
end
blnk=' ';
fmt1a='%s';
fmt1e='%s\n';
fmt1b='%6.3f%6.3f%6.3f%6.3f%6.3f%6.3f';
fmt1c='%6.3f%6.3f%6.3f%6.3f%6.3f%6.3f\n';
fmt1d=[fmt1b fmt1c];
for n=1:nreg;
	x=RS(n,:);
	fprintf(fid3,'#%3.0f',n);
	fprintf(fid3,fmt1d,x);
end
fclose(fid3);


% Station sample size (#yr avail in reference period) on which
% unadjusted means and std devs are based
file3=['yrsref' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
blnk=' ';
fmt1a='%s ';
fmt1e='%s\n';
fmt1b='%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f';
fmt1c='%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n';
fmt1d=[fmt1b fmt1c];
% cull non-NaN rows
X=USIZE;
L1=isnan(X);
L2=(~all(L1'))';
X=X(L2,:);
mX=size(X,1);
BB1=B1(L2,:);
for n=1:mX;
	x=X(n,:);
	fprintf(fid3,fmt1a,BB1(n,1:5));
	fprintf(fid3,fmt1d,x);
end
fclose(fid3);


% The ratios of monthly means used to adjust station means 
file3=['ratmn' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
fmt1b='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f';
fmt1c='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n';
fmt1d=[fmt1b fmt1c];
% cull non-NaN rows
X=RATMN;
L1=isnan(X);
L2=(~all(L1'))';
X=X(L2,:);
mX=size(X,1);
BB1=B1(L2,:);
for n=1:mX;
	x=X(n,:);
	fprintf(fid3,fmt1a,BB1(n,1:5));
	fprintf(fid3,fmt1d,x);
end
fclose(fid3);


% The ratios of monthly std devs used to adjust station std devs 
file3=['ratsd' filecode '.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
fmt1b='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f';
fmt1c='%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n';
fmt1d=[fmt1b fmt1c];
% cull non-NaN rows
X=RATSD;
L1=isnan(X);
L2=(~all(L1'))';
X=X(L2,:);
mX=size(X,1);
BB1=B1(L2,:);
for n=1:mX;
	x=X(n,:);
	fprintf(fid3,fmt1a,BB1(n,1:5));
	fprintf(fid3,fmt1d,x);
end
fclose (fid3);


% Number of good stations matrix, all regions combined
file3=['nstns' filecode '.dat'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
blnk=' ';
fmt1a='%s ';
fmt1e='%s\n';
fmt1b='%4.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f ';
fmt1c='%3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n';
fmt1d=[fmt1b fmt1c];
X=[yr3 NSTNS];
fprintf(fid3,fmt1d,X');
fclose(fid3);


% Individual regional files of number of good stations over time
% Note that the counts in NSTNS and NS apply to stations as crosslisted
% in IR1, and the the regional climate series might cover a longer
% period of data.  This is because the regional climate series can
% be built with distance weighting from stations "outside" the region.
% For example, the "south" san pedro basin might use Tucson for 
% regional climate up until Benson comes in the early 1880's, while
% NS for south san pedro would say no stations before Benson. The regional
% NS do give useful info on reliability of regional info -- saying,
% for example, that the pre-1880s san pedro regional series is really
% interpolated from long distance.

for n=1:nreg;
	% column pointer to NS
	jgo=1+(n-1)*12;
	jsp=jgo+11;
	jcols=jgo:jsp;
	% get regional count
	X=NS(:,jcols);
	yr3a=yr3;
	% trim all years with all elements 0
	L1=X==0;
	L2=(all(L1'))';
	if sum(L2)>0;
		X(L2,:)=[];
		yr3a(L2,:)=[];
	end
	d1=diff(yr3a);
	if ~all(d1==1);
		disp(['region ' int2str(n)]);
		error('yr3a for above region not continuous');
	end
	% build output file
	flout=['nstns' filecode int2str(n) '.dat'];
	pf3=[path5 flout];
	fid3=fopen(pf3,'w');
	X=[yr3a X];
	fprintf(fid3,fmt1d,X'); % note that use fmt1d build in prev section
	fclose (fid3);
end

fclose all;


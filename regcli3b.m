function regcli3b(flsin,miscin,Sinf,pntnum)
%
% Last revised 7-24-98
%
% Standardized anomaly method for regional climate series -- point-by-point
% version.  Adapted from regcli2.m  Regcli3b.m interpolates station PPT or TMP 
% monthly data to gridpoints, or to any network of lon-lat points.  Big
% differences from regcli2.m:
%
%	-regcli3bb has no user prompted info; info is passed as arguments
%	-regcli3b does not implement adjustment of reference period means or standard
%		devs. In terminology of regcli2.m, a 'pass-1' analysis is all that is
%		done
%	-regcli3b is intended to be called for individual gridpoints (or, say, tree
%		tree-ring sites) for which interpolated series are desired.  For a 
%		western North America application, calling function would call regcli3b.m
%		repeatedly to build up interpolated series
%	-distance computations are done in calling function regcli3a
%  -regcli3b has a modification (7-24-98) that substitutes zero as the interpolated
%     pcp for any month/year in which the pcp was zero at all predictor stations.
%     The zero overrides the interpolated value.  All sections with this
%     modification are flagged with "Revision 7-24-98"
%
% Application.  Wanted tree-site-centered TMP series at 700+ western NA tree
% sites, from network of 1002 GHCN v2 TMP stations. Direct use of regcli2.m 
% for this problem would be difficult because of huge array sizes. For example,
% 1002 x 700 just to store the station-to-tree distances, and 
% 100 yr x 700 stations x 12 months to store all stations' monthly data.
% Decided instead to write regcli3b.m to treat problem tree-site by tree-site
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
%     {8} kmassive (1 x ?)s  "Yes", indicating a massive generic run to get
%         all the urban heat island adusted series in files TD1.mat, TD2.mat, etc
%         or "No", indicating a more general type of run to maybe interpolate
%         pcp or tmp to points
%     {9} preftype (1 x 1)s  'p' for pcp, 't' for tmp.  If kmassive "Yes", preftype
%         is empty.  Note that if preftype is p, outfiles will be named
%         pntp1, pntp2, etc.  If preftype is t, names are pntt1, pntt2, ...
%	Sinf {} string matrices of identifying information
%		{1} (1 x 70)s for the target point
%		{2} (? x 70)s for the climate stations represented by the .mat files
%
%	pntnum(1 x 1)i  point number -- numeric identifier used to build prefixes
%			of output files
%
%**************************** OUTPUT FILES ***********************************
%
% A .mat file holding the target point monthly climate series in matrix Z
% The name of this .mat file is built from the input variable pntnum, as
% ['P' int2str(pntnum)] for pcp, and ['T' int2str(pntnum) for tmp
%
% A .mat file holding the sample size (number of stations) going into the
% target point interpolation.  Same size as the .mat file of monthly data
%
% An ascii information file for the point, with same prefix as the .mat file,
% but with .txt extension.  The first line of this file contains information
% on the target point itself: 
%		Sequence number of point
%		Prefix of .crn file 
%		tree-ring site name (if a tree-ring site)
%		species
%		lon and lat
%		elevation in m
%		first and last year with any monthly data in the interpolated series
%
% Lines 2-6 contain information on the climate stations used
% to form the target-point reference-period means and standard deviations.
% These are usually the 5 nearest stations, or sometimes fewer than 5 depending
% on kopts(3). Distance alone is not the only criterion.  The station record
% must also have the minimum required percentage of valid data for each month
% in the reference period. Information on lines 2-6 is:
%		seq number (nearest as #1); 
%		Station .mat file name prefix
%		station name
%		country
%		lon/lat, 
%		elevation, 
%		time covrage (first,last yrs with any monthly data at station)
%		distance to target point
%		weight for the target-point mean and std dev
%
% Example:
% 
%	1-az545 		Tsegi Point PIED  	-110.32	34.56  1243m 1904 1995
%	1-D67			Holbrook 	USA    	-109.23 35.64   934m 1912 1996 32km  0.123
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
% wgtdist2.m inverse distance weighting
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


%***************   Unload cell information

CG = miscin{1};
C2 = miscin{2};
A = miscin{3};
yrinf=miscin{4};
pdref=yrinf(1,:);
pdstore=yrinf(2,:);
kopts=miscin{5};
kmode=kopts(1);

dcrit = miscin{6}(1);
dsrch = miscin{6}(2);
nmax = miscin{6}(3);

kmassive = miscin{8};
preftype=miscin{9};


%********************************************************************
% Get the string matrix of filename of monthly climate stations, 
% with path attached
B1=flsin;
ns1 = size(B1,1); % number of input .mat climate files


%*********** QC the input args
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
end
if nC2~=2;
	error('C2 should be 2-column');
end
if  any(any (C2(:,1)>0));
	error('C2 long should be negative')
end
if  any(any(C2(:,2)>90));
	error('C2(:,2) should be all 90 deg or less (latitude)');
end
clear mC2 nC2


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
%disp('READING MONTHLY DATA FILES AND STORING DATA');
%disp(' ');
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
	eval(['load '  fln]);
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

%disp('WEIGHTING STDZD ANOMALIES TO TARGET POINT');

% Allocate
P=a(ones(nrows,1),ones(12,1)); % to hold target point data in clim units
Q=a(ones(nrows,1),ones(12,1)); % to hold target point data in stdzd units
NS=a(ones(nrows,1),ones(12,1)); % to hold sample size (number of stations
% weighting point) in each month/year

%revision 7-24-98
% Make logical pointer to years/months with precipitation zero at all 
% available stations used as predictors.  Will use this information to
% substitute zero for the interpolated precip for those years/months.
% Lallzero will not be used if the climate variable is not pcp
Lallzero=logical(zeros(nrows,12));
% revision end
   
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
      
	%disp(['   Month ' int2str(kmon)]);
	jp2=JP2(kmon); % scalar index to col of P,Q for this month

	% Pull the subset of station stdzd anomalies, all yr in E 
	F1=E(:,J2(:,kmon));
   
      
	% Pull subset of rows of F1 for which not all values are NaN
	L1=~isnan(F1);
	L1a=(any(L1'))';
	F2=F1(L1a,:);
   
   % Revision 7-24-08
   % Pull subset of station data in original units (all years of this data are in D)
   % Need these data to check for years in which precip is zero at the only 
   % available predictor stations.
   F3 = D(:,J2(:,kmon)); % all years
   F4 = F3(L1a,:); % just the years for which not all the predictors are NaN.
   % F4 is the orig data corresponding to the stdzd anomalies in F2
   % Revision end
   
   
	% Make logical pointer to F2 with data
   L1=~isnan(F2);
   
   % Revision 7-24-98
   % Make pointer to rows (years/this month) where data values for 
   % all stations either zero or NaN.  Note that already have scrrened so that
   % no row of F2 and F4 is all NaN
   Ltemp = F4==0 | isnan(F4);
   LZall=(all(Ltemp'))'; % col vector. Is 1 if the only available predictor data
        % precip in the month/year is zero. Is 0 otherwise
   % Revision end
   
   
	% Compute matrix of weights for this month of year
	kk=[kopts(3) kopts(4)];
	[WW,NW,JW]=wgtdist2(A,L1,dcrit,dsrch,nmax,kk);
      
	jpthis=jp2; % col in P, Q for this month		

	colsWW=JW(1):JW(2);  % columns of W with this points wts
	W1=WW(:,colsWW);
		
	% Weight the stations
	G1=F2 .* W1;
      
	% Sum over cols to get cv of weighted stdzd anoms this point/month
	g1 = (sumnan(G1'))';

% Store weighted anomalies and sample size
  Q(L1a,jpthis)=g1;
  NS(L1a,jpthis)=NW;
  
  % Revision 7-24-98
  % Store information on whether all predictor values zero (used for pcp only)
  Lallzero (L1a,jpthis)=LZall;
  % Revision end
  
  
end



%************** GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS
%
%disp('COMPUTING GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS');

% Allocate
PM=a(:,ones(12,1)); % to hold means
PS=a(:,ones(12,1)); % to hold std devs
W5=a(:,ones(ns1,1)); % to hold weights on stations

TMN=UMN;  % means
TSD=USD; % standard devs


% Mean
L1=~isnan(TMN');
[WWM,NWM,JWM]=wgtdist2(A,L1,dcrit,dsrch,nmax,kk);
% Note that all rows of WWM are identical. Only need first
for n=1:mG;
	jcols=JWM(n,1):JWM(n,2);
	W5(n,:)=WWM(1,jcols);
end
L1=isnan(W5);

% Store pointer to selected stations for weighting means, and the weights
Lstns = ~L1;
w5=W5(Lstns);

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


%------  Next two operations aimed at pcp data.  Possible that have 1 or more
% predictor stations in a year/month, but that pcp zero at all those stations. 
% Blind application of stdzd anomaly method would typically treat zero as some
% number of standard deviations below the ref-period mean for the station. 
% Weighting the stdzd anoms over stations would then give perhaps some small
% negative number as the mean standardized anomaly for the month/year.  The
% building back in the regional standard deviation and mean could make the 
% final interpolated value, in pcp units, either (1) negative, or (2) some small
% positive value.  But negative pcp is illogical.  And it seems reasonable that
% if pcp is zero at all existing stations, best to give zero as the estimate
% for the gridded value. Thus the next two operations.

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


%Revision 7-24-98
% If data pcp, replace interpolated value with zero if all predictor stations
% had zero pcp in year/month
if kopts(2)==1; 
   P(Lallzero)=0;
elseif kopts(2)==2; % "no-zero-substitute mode", appropriate for tmp 
   % No action
else
   error('Invalid setting for kopts(2)');
end
% Revision end

%------ end of the two special operations for pcp or pcp-like variables



%******************** *****************************
%disp('BUILDING OUTPUT FILES')
 
%****** OUTPUT THE RESULTS **************************
%
% 
% RECALL THAT
%
% B1 station filename list
% C2 long lats of stations
% CG long lats of point
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



% Reduce row size of P, Q, NS by lopping off any trailing or leading
% years that are all NaN
L1=isnan(Q);
L2=(all(L1'))';
irow=(1:length(L2))';
irow(L2)=NaN;
irow=trailnan(irow);
irow=flipud(irow);
irow=trailnan(irow);
irow=flipud(irow);
irowgo = irow(1); irowsp = irow(length(irow)); % start, end indices of series
yr4 = (yr2(irowgo):yr2(irowsp))';
nyrs4=length(yr4);
P1 = P((irowgo:irowsp),:);
Q1 = Q((irowgo:irowsp),:);
NS1 = NS((irowgo:irowsp),:);

% Store first and last year for which any data exists for the point, and
% for the continuous period with no missing data
yrgo1=yr4(1);  yrsp1=yr4(nyrs4);
[yrgo2,yrsp2]=nonan1(P1,yr4);



% Build and save individual  .mat  files of data and sample sizes
pathout=miscin{7}(2);
pathout=pathout{1};
X = [yr4 P1];
N = [yr4 NS1];

if strcmp(kmassive,'Yes');
   flout1 = [pathout 'TD' int2str(pntnum)];
   eval(['save ' flout1 ' X;']);
   flout2 = [pathout 'N' int2str(pntnum)];
   eval(['save ' flout2 ' N;']);
else
   flout1 = [pathout 'pnt' preftype int2str(pntnum)];
   eval(['save ' flout1 ' X;']);
   flout2 = [pathout 'N' preftype  int2str(pntnum)];
   eval(['save ' flout2 ' N;']);
end


% Build and save ascii file of point information; want like
%	1-az545 		Tsegi Point PIED  	-110.32	34.56  1243m 1904 1995 1915 1990
%	   1-D67			Holbrook 	USA    	-109.23 35.64   934m 1912 1996 32km  0.123
%

if strcmp('kmassive','Yes');
   pftxt = [pathout 'info.txt'];
else
   pftxt = [pathout 'info' preftype '.txt'];
end

fid8 = fopen(pftxt,'a');

distall=A(Lstns);
stnall = Sinf{2}(Lstns,:);
% Find max length of any line of stnall
maxlen=0;
for n =1:length(w5);
   c=deblank(stnall(n,:));
   maxlen=max(maxlen,length(c));
end
eval(['fmt1= ''%' int2str(maxlen) 's %4.0f  %6.5f'';']);  
   
str1=deblank(Sinf{1});
str1=sprintf('%s %4.0f-%4.0f   %4.0f-%4.0f',str1,yrgo1,yrsp1,yrgo2,yrsp2);
fprintf(fid8,'%s\n',str1);
for n =1:length(w5); % loop over selected sites
   dist=distall(n); % point to station distance (km)
   w5this =w5(n); % weight on station
   str2 = deblank(stnall(n,:));
   str2 = sprintf(fmt1,str2,dist,w5this);
   fprintf(fid8,'%s %s\n',blanks(1),str2);
end
fprintf(fid8,'%s\n',blanks(5));

fclose all;


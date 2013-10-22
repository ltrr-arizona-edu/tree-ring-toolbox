function regcli4
% regcli4: regional P or T from station data by stdzd normal method-- SRER version
% regcli4;
% Last revised 5-22-99
% 
%*** STEPS
%
% 1) Load info from fls?.mat:
%  path3 (1 x ?)s path to input .mat monthly climate files
%  path5 (1 x ?)s path to where to put output
%  names1 {? x 1}cell of string name codes for the stations
%  C2 (? x 2)r  long and lat for stations in names1; decimal deg
%  IM (? x 1)i row index to names1 tells master stations
%  iuse (? x 1)i index to rows of names1 telling which stations eligible for 
%     use in the interpolation (iuse created by threegd.m)
%
% 2) Load info from misc?.mat:
%  CG (? x 2)r lon, lat of gridpoints
%  IR3(? x 1)i row index to CG tells gridpoints to use in making regional series
%  dcrit (1 x 1)r  minimum distance used in weighting stations (km) 
%  dsrch (1 x 1)r  search radius for looking for stations around gridpoints
%    (ignored if kopts(3)==1)
%  nmax (1 x 1)i maximum number of stations to use in weighting any gridpoint\
%    (see kopts(3))
%  pdref (1 x 2)i  first, last year of reference period for adjusting means
%  pdstore (1 x 2)i first last year of data-storage period -- must cover all
%      years for which output data to be stored
%  B (2 x 12)r  regression constant (row 1) and coef for std dev vs mean, each month
%  dataZ (1 x 1)s the variable (e.g., X) storing monthly climate data in the monthly .mat files
%  kopts(1 x 4)i options
%     kopts(1) adjustment method 
%       ==1 ratio method for mean, line of std dev vs adjusted mean for std dev
%       ==2 no adjustment; use observed station means and standard deviations
%           (note that need 2 observations to get a mean and std deviation)
%       ==3 adjust mean as in (1), but use observed standard deviation
%     kopts(2) handling of negative interpolated values
%       ==1 assign as zero (relevant to pcp)
%       ==2 accept negative  (relevant to temperature)
%     kopts(3) search rule
%       ==1 disregard search radius and use the nearest nmax stations for weighting;
%           if fewer than nmax stations available in any given month, use that number
%       ==2 require that stations be in search radius, and use the nearest 
%           up-to-nmax of those in the weighting. May need trial and error runs to
%           make sure search radius large enought
%     kopts(4) interpolation method
%       ==1 inverse distance (only implemented method so far)
%    

% Template was regcli2.  Wrote regcli4.m specifically to handle McClaran's 
% Santa Rita Expt Range precip data.  Differences from regcli2.m:
%
% 1) straight line reltn std dev vs mean used to estimate monthly standard
%    deviations for some stations.
% 2) mean-ratio OK for estimating adjusted monthly mean even if only one
%    observed year of data for a given month.  Method is to estimate the
%    adjusted mean first, then to use the straight-line fit vs the
%    adjusted mean (from long stations) to get the estimated standard
%    deviation.  Standardized departure is (O-M)/s, where is O is the 
%    monthly P, M is the adjusted mean for that month, and s is the
%    standard deviation (computed if from one of the long-term stations,
%    and estimated if from another of the stations).
% 3) User prompted to point to file holding .mat station filenames and also:
%    - pointer to the long-term stations
%    - coef matrix allowing estimation of monthly std dev for any station/month,
%      given the adjusted mean.
%
% Weighting implemented so far is  inverse distance
% 
%*** INPUT FILES *********************************
% 
% 1. fls*.mat input .mat file with the following stored variables
%     path3 (1 x ?)s path to monthly climate .mat files (e.g., 'c:\projs\ai6\stnpcp\')
%     path5 (1 x ?)s path to output files to be created (e.g., 'c:\projs\ai6\outfls\')
%     names1 (? x 1)cell:  idcodes of stations (all 75) 
% 2. misc*.mat -- miscellaneous program control
%		CG [mG x 2]r long-lat file in mapping units for the mG gridpoints.
%		C2 [ns1 x 2]r long-lat file in mapping units for the ns1 stations
%		IR2 [? x 1]i index to rows of CG telling which gridpoints are
%			to be used in computing the regional series
%		IM (? x 1)i  index to rows of C2 telling which stations are to
%			be used as master stations
%		pdref (1 x 2)i reference period (start, end years) for adjusting
%			means
%		pdstore (1 x 2)i  start, end years for storage matrix that must
%			cover at least the period with any valid data at any station.
%
%			dcrit -- critical threshold distance to avoid over-weighting site in
%				distance weighting (km)
%			dsrch -- search distance for weighting (km)
%			nmax --  max number of stations to weight to a gridpoint
%

% 
%**************************** OUTPUT FILES ***********************************
%
% The "?" below refers to a screen prompted letter code for the run
% The i refers to region i
% The j refers to gridpoint j
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
% 
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
% Reference: Jones & Hulme (1996).  The reference does not cover the case
% of adjusting the means and standard deviations of station data to
% account for anomalous climatological reference periods.
%
% Method.   Monthly station time series are first converted to standardized
% anomaly series using the mean and standard deviation for that station/month.
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
%% The standardized anomaly series are then weighted over stations around a 
% gridpoint  (subject to constraint on search radius) to produce a gridpoint
% average standardized anomaly monthly time series.
%
% The gridpoint standardized anomalies are converted to original units
% (e.g., in for ppt) using a weighted average of reference-period means
% and standard deviations for stations in the search radius, dsrch, around
% the gridpoint.
%
% The regional climate series is computed by averaging the gridpoint series
% over the gridpoints specified by IR2.  Note that the averaging is done 
% over gridpoints series after they have been converted to original units.
% Gridpoints are weighted proportional to cosine of the latitude in 
% producing the regional series to adjust for convergence of the meridians.
%
% "Gridpoints" might be long-lat points in application to produce regional
% climate series.  Gridpoints might also be tree locations in an application
% to produce tree-location-weighted climate series.  In this case, the 
% gridpoints are irregularly spaced.
%
%
%

%-- Prompt for just the facts, or expanded output mode
% kmode 1 is full-output model. kmode 2 is just the facts
kfull=questdlg('Expanded-output mode?');
switch kfull;
case 'Yes';
   kmode=1;
otherwise;
   kmode=2;
end;
clear kfull;
% kmode quick vs diagnostic mode

%--- Load path and station info
[file1,path1]=uigetfile('fls*.mat','Infile with path and station info');
pf1=[path1 file1];
eval(['load ' pf1  ' names1 path3 path5 C2 IM  iuse dataZ;']);
% Check that variables exist
Ltemp1=[exist('names1') exist('path3') exist('path5') exist('C2') exist('IM') exist('dataZ')];
Ltemp1a=[exist('iuse')];
Ltemp1=[Ltemp1 Ltemp1a];
Ltemp=all(Ltemp1==1);
if ~Ltemp
	clc
	disp(['The following variables should be in ' pf1]);
	disp('names1 path3 path5 C2 IM iuse dataZ');
	error('All the above were not in the file');
end
clear Ltemp Ltemp1;
% file1, path1, pf1 -- path and station info

%---- Load other input information from misc?.mat
[file2,path2]=uigetfile('misc*.mat','Infile with miscellaneous info');
pf2=[path2 file2];
eval(['load ' pf2  ' B CG  IR2 dcrit dsrch nmax kopts pdref pdstore;']);
% Check that variables exist
Ltemp1=[exist('B') exist('CG') exist('IR2') exist('dcrit') exist('dsrch')];
Ltemp2=[exist('kopts') exist('pdref') exist('pdstore') ];
Ltemp3=[Ltemp1 Ltemp2];
Ltemp=all(Ltemp3==1);
if ~Ltemp
	clc
	disp(['The following variables should be in ' pf2]);
	disp('B CG IR2 dcrit dsrch nmax kopts pdref pdstore');
	error('All the above were not in the file');
end
clear Ltemp Ltemp1 Ltemp2 Ltemp3;
% file2,path2, pf2    path & file for misc?.mat input


%------ QC of Input

if ~iscell(names1);
   error('names1 should be cell');
else;
   ns1 = size(names1,1); % ns1 is number of stations in names1
end

if size(C2,1)~=ns1;
   error('row size of C2 must match that of names1');
end
Ltemp=[size(C2,2) size(CG,2)];
if ~all(Ltemp==2);
   error('C2 and CG must have column size of 2');
end;
if any(C2(:,1)>0) | any(CG(:,1)>0);
   error('C2 and CG longitudes must be negative');
end
if any(C2(:,2)>90) | any(CG(:,2)>90);
   error('C2 and CG latitudes greater than 90 deg');
end
[mG,nG]= size(CG); % number of gridpoints is mG
clear Ltemp;


% B has the regression coefs of std dev vs mean for each month
[mB,nB]=size(B);
if kopts(1)==1;
   if (mB~=2) | (nB~=12);
      error('B must be 2 x 12');
   end
else;
   if ~isempty(B);
      error('B should be empty if kopts(1)~=1');
   end
end


% IM says which stations are to be used in interpolating the
% regional master series 
[nmast,ntemp]=size(IM); % nmast is the number of stations in master
if(ntemp~=1);
   error('IM should be col vector');
end
if any(IM>ns1);
	error('An element of IM larger than highest station seq no.');
end
clear ntemp;

% iuse tells which stations in names1 are eligible for use in the analysis
nuse = length(iuse);  % number of stations to use
% Check that all master stations are included in the set marked by iuse
for n = 1:nmast;
   itemp = IM(n);
   nmtemp = names1{itemp};
   Ltemp=itemp==iuse;
   if sum(Ltemp)~=1;
      error(['master station ' nmtemp ' not in iuse list']);
   end
end;
clear itemp nmtemp Ltemp;


% IR2 says which gridpoints are to be used in computing the regional series
[npts1,nIR2]=size(IR2);
if npts1>mG | nIR2~=1;
   error('IR2 must be col vector with no more elements than CG has rows');
end
if any(IR2>mG)
	error('An element of IR2 larger than highest grdpnt seq no.');
end
% npts1 is number of gridpoints to be used in regional average

% pdref holds start, end years for reference period
if pdref(1)>pdref(2); 
	error('Reference period ends before it starts');
end
yr1=(pdref(1):pdref(2))'; % year vector for reference period
nref=length(yr1); % # yrs in ref pd

% pdstore holds start, end years of a matrix to store the monthly
% data for all stations. 
if pdstore(1)>pdstore(2)
	error('end year of pdstore precedes start year');
end

% kopts holds program options
if kopts(1)<1 | kopts(1)>3 
	error('Illegal setting of kopts(1)')
end
if kopts(2)<1 | kopts(2)>2;
	error('Illegal setting of kopts(2)')
end
if kopts(3)<1 | kopts(3)>2;
	error('Illegal setting of kopts(3)')
end
if kopts(4)<1 | kopts(4)>1;
	error('Illegal setting of kopts(4)')
end



%-- RE-EXPRESS INDEX IM IN TERMS OF REDUCED SET OF nuse STATIONS; ALSO
%  REDUCE OTHER VARIABLES TO COVER ONLY THOSE nuse STATIONS
Lt1 = zeros(ns1,1);
Lt1 (iuse)=1;
I1 = cumsum(Lt1);
IM1 = I1(IM);
names1=names1(iuse,:);
IM=IM1;
C2=C2(iuse,:); % long / lat of stations
clear IM1 I1 Lt1;



% Optional feedback on station composition in regions
if kmode==1; % diagnostics mode
	clc
	disp('Listing of all monthly climate stations')
	
	disp(names1);
	
	disp('Press any key to continue')
   pause
   clc
	
	disp('Data:');
   disp([int2str(ns1) ' stations in names1']);
   disp([int2str(nuse) ' stations to be used in full interpolation ']);
   disp([int2str(nmast) ' stations in master']);
   disp([int2str(mG) ' gridpoints'])
   disp([int2str(npts1) ' gridpoints used in regional average']);
   
   disp('Reference period for master series')
	disp([pdref(1) pdref(2)]);
	
	disp('Period for huge matrix');
	disp([pdstore(1) pdstore(2)]);
	
	disp('Settings for call to wgtdista.m')
	disp([' dcrit = ' num2str(dcrit) ' km'])
	disp([' dsrch= ' num2str(dsrch) ' km'])
   disp([' nmax= ' int2str(nmax) ' stations'])
   disp([' Option for method: ' int2str(kopts(1))]);
   disp([' Option for negatives: ' int2str(kopts(2))]);
	disp([' Option for search : ' int2str(kopts(3))])
	disp([' Option for weighting method: ' int2str(kopts(4))])
		
	disp('Press any key to continue')
   pause
   clc
end


%--COMPUTE MATRICES OF DISTANCE FROM GRIDPOINTS TO STATIONS

% get index to stations for master
im=IM;
Dist1 = gcdist(CG,C2); % gridpoints to all stations
Dist2 =  Dist1(:,im); % gridpoints to master stations


% Allocate for adjusted and unadjusted means and std devs, and for adjustment ratio
AMN=repmat(NaN,nuse,12);
ASD=repmat(NaN,nuse,12);
UMN=repmat(NaN,nuse,12);
USD=repmat(NaN,nuse,12);
RATMN = repmat(NaN,nuse,12);
USIZE=repmat(NaN,nuse,12);


%**************** READ IN MONTHLY STATION CLIMATE SERIES
%
disp('READING MONTHLY DATA FILES AND STORING DATA');
disp(' ');
% Store the input monthly ppt series for all specified climate stations for
% the region in a single matrix

% Initialize matrices 

yr2=(pdstore(1):pdstore(2))'; % years vector for storage matrix
nrows = length(yr2); % number of rows (years) in storage matrix
ncols = nuse*12; % % number of cols in storage matrix
D=repmat(NaN,nrows,ncols); % raw monthly data
E=repmat(NaN,nrows,ncols);  % monthly data, unadjusted anoms
F=repmat(NaN,nrows,ncols);  % monthly data, adjusted anoms
T=repmat(NaN,nuse,2);  % start, end years of each station's record
I=repmat(NaN,nuse,2);  % corresp row pointer to D


% What cols in D and E will the station monthly data go into?
j1=(1:12:(ncols-11)); % start col for stn1, 2,etc
j2=j1+11; % end col
J1=[j1' j2'];
% First row of J1 will give first and last storage row for stn 1, second row
%		stn 2 etc

% What cols in D and E will hold all Jan data, Feb data, etc?
j1=(1:12);
j1=j1(ones(nuse,1),:);
j2=(1:nuse)';
j2=12*(j2-1);
j2=j2(:,ones(12,1));
J2=j1+j2;
% First col of J2 will point to cols of D for Jan, secnd for Feb, etc

% Loop over stations: get first and last year,
% 
for n=1:nuse;
	clear Z;
	fln=names1{n};
   eval(['load ' path3 fln]);
   eval(['Z = ' dataZ ';']);
   eval(['clear ' dataZ ';']);
   
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


%****************************************************************
% Compute the master series for the reference period, and the 
% reference-period monthly means and standard deviations for the
% master series. 

% Rename number of stations in regional series and in master
num2=nuse; % # stns in region
num3=nmast; % # stns in master
B2=names1(im,:); % cell of station ids for master
B3=names1; % stations in region 
if kmode==1;
   disp('Series specified for master')
   disp(B2);
   disp('Press any key to continue')
   pause	
end


%----- MASTER STATION ANALYSIS-- ALLOCATE
% Recall that num3==nmast is the number of master stations
numcol=num3*12; % # cols in master storage mtx
M=repmat(NaN,nref,numcol); %  ref pd raw data, years by months
LL=repmat(NaN,nref,numcol); %  ref pd missing value pointer
V=repmat(NaN,nref,numcol); % ref period station stdzd anomalies
VM=repmat(NaN,nref,12); % mean stdzd anomaly series for the master
YM=repmat(NaN,nref,12); % VM converted to original climate units
Mrefm=repmat(NaN,num3,12); % ref period long-term monthly means
Mrefs=repmat(NaN,num3,12); % ref period long term std devs
gm=repmat(NaN,1,12); % global mean monthly ppt (average over stations)
gs=repmat(NaN,1,12); % global mean std dev

%--- MASTER STATIONS -- COMPUTE STATION STDZD ANOMS AND STORE MONTHLY DATA
% Get the monthly data for the master stations and store in matrix M
Lmast = yr2>=pdref(1) & yr2<=pdref(2); % ref pd rows in the mtx
for nn = 1:num3; % loop over master stations
   j1=J1(im(nn),:); % start,end cols for this station in the huge mtx
   kcols=j1(1):j1(2);
   Z=D(Lmast,kcols); % the monthly climate data for ref period
   % Compute the long-term monthly means and std deviations
   zmn = nanmean(Z);
   zst = nanstd(Z);
   
   % Compute standardized anomalies
   ZMN = repmat(zmn,nref,1);
   ZST = repmat(zst,nref,1);
   V1 = (Z - ZMN) ./ ZST;
   
   % Compute storage columns for this data in M
   jgo = (nn-1)*12 +1;
   jsp = jgo+11;
   jcols=[jgo:jsp];
      
   % Store results
   M(:,jcols)=Z; % monthly data
   V(:,jcols)=V1; % standardized anomalies
   Mrefm(nn,:)=zmn; % ref-pd monthly means
   Mrefs(nn,:)=zst; % ref-pd monthly std devs
   
   % Compute global mean and std dev
   gm=mean(Mrefm);
   gs=mean(Mrefs);
end; % for nn=1:num3
% M monthly climatic data (yrs x (mos x stns)) for master stations
% V standardized anomalies corresponding to M
% Mrefm, Mrefs  ref-pd long-term (unadjusted) means and std devs for master stations



%--- MASTER:  CHECK THAT EVERY MONTH/YEAR HAS ADEQUATE COVERAGE IN MASTER
% 
% Require at least that each month/year be represented by at least 2 stations.

% What cols in M and V hold all Jan data, Feb data, etc?
j1=(1:12); % a row vector
j1=repmat(j1,num3,1); % duplicate to row-size equal to number of master stations
j2=(1:num3)'; % col vector  [1 2 3  ... num3]'
j2=12*(j2-1); % col vector [0 12 24 36 ... (num3-1)*12]'
j2=j2(:,ones(12,1)); % col-dupe to 12 identical cols. Size now num3 x 12
JM2=j1+j2; % pointer to columns of M and V
clear j1 j2;
% First col of JM2 will point to cols of M with Jan data for all
% num3 master stations.  Second col to cols of M with Feb data, etc

%--- Loop over months
for n =1:12;
   j1 = JM2(:,n); % col reference to M for month n
   MM = M(:,j1); % subset of ref-pd master station monthly data for all stations, this month
   L1 = ~isnan(MM);
   L1sum = sum(L1,2); % col vector of number of valid values in each year for this month
   if any(L1sum<2);
      error(['Fewer than two data values in at least one year in month ' int2str(n)]);
   end;
end;
clear j1 MML1 L1sum


%-- MASTER: INTERPOLATE MASTER STATION DATA TO GRIDPOINTS

% Allocate
nP = mG*12; % #cols in huge gridpoint tsm
GPTm=repmat(NaN,mG,12); % gridpoint long-term means interpolated from master stations
GPTs=repmat(NaN,mG,12); % gridpoint long-term std devs interp from master stns
GPT1a=repmat(NaN,nref,nP); % gridpoint interp monthly stdzd anoms
GPT1 =repmat(NaN,nref,nP); % gridpoint interp monthly values

%---- Build pointers to columns of GPT1a and GPT1 for Jan, Feb, etc.
j1=(1:12:(nP-11)); % start col for point 1, 2,etc
j2=j1+11; % end col
JP1=[j1' j2'];
% First row of JP1 will give first and last storage row for 
% point 1, second row for point 2 etc.  Points numbered as in 
% input coordinate mtx CG

% What cols in GPT1a and GPT1 will hold all Jan data, Feb data, etc?
j1=(1:12);
j1=j1(ones(mG,1),:);
j2=(1:mG)';
j2=12*(j2-1);
j2=j2(:,ones(12,1));
JP2=j1+j2;
% First col of JP2 will point to cols of GPT1a or GPT1 with Jan data for all
% mG gridpoints.  Second col to cols of GPT1a or GPT1 with Feb data, etc


% Loop over the 12 months of the year
for kmon=1:12;
   disp(['   Month ' int2str(kmon)]);
   
   jp2=JP2(:,kmon); %cv of target cols indices to GPT1a, GPT1 for this month
	jm2=JM2(:,kmon); %cv of source cols indices to stdz anomalies V for this month

	% Pull the subset of station stdzd anomalies, ref pd years
	F1=V(:,jm2);
	
	% Pull subset of rows of F1 for which not all values are NaN
	L1=~isnan(F1);
	L1a=(any(L1'))';
	F2=F1(L1a,:);

	% Make logical pointer to F2 with data
	L1=~isnan(F2);

	% Compute matrix of weights for this month of year
	kk=[kopts(3) kopts(4)];
	[WW,NW,JW]=wgtdist1(Dist2,L1,dcrit,dsrch,nmax,kk);

	% Loop over the gridpoints
	for kgrd=1:mG;
		if kmode==1;
			disp(['GRIDPOINT # ' int2str(kgrd)]);
		end

		jpthis=jp2(kgrd); % col in GPT1a, GPT1 for this month/gridpoint		

		colsWW=JW(kgrd,1):JW(kgrd,2);  % columns of W with this points wts
		W1=WW(:,colsWW);
		
		% Weight the stations
		G1=F2 .* W1;
      
		% Sum over cols to get cv of weighted stdzd anoms this point/month
		g1 = (sumnan(G1'))';

		% Store weighted anomalies
		GPT1a(L1a,jpthis)=g1;
	end
end


%--  WEIGHT LONG-TERM MEANS AND STD DEVS OF MASTER STATIONS TO GRIDPOINTS

disp('COMPUTING MASTER-TO-GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS');

% Note that GPTm and GPTs, both nmastx12, will hold these results
% Note that Mrefm and Mrefs, both nmast x 12, hold the long-term master stn means and stds

W5M=repmat(NaN,mG,nmast); % will hold weights on master stations


% Mean
L1=~isnan(Mrefm');
[WWM,NWM,JWM]=wgtdist1(Dist2,L1,dcrit,dsrch,nmax,kk);
% Note that all rows of WWM are identical. Only need first
for n=1:mG;
	jcols=JWM(n,1):JWM(n,2);
	W5M(n,:)=WWM(1,jcols);
end
L1=isnan(W5M);
sum5=sum(sum(L1));
% Matrix mult with NaNs gives Nans, so change NaNs to zero
if sum5>0;
	W5M(L1)=zeros(sum5,1);
end

% Make temporary matrices Hm and Hs by replacing NaNs with zeros
% in TMN, TSD before matrix mult.  Otherwise, NaNs in matrix mult, which
% means all products NaN
Hm=Mrefm;
Hs=Mrefs;
L1=isnan(Hm);
sum5=sum(sum(L1));
if sum5>0;
	Hm(L1)=zeros(sum5,1);
end
L1=isnan(Mrefs);
sum5=sum(sum(L1));
if sum5>0;
	Hs(L1)=zeros(sum5,1);
end

GPTm=W5M*Hm;
GPTs=W5M*Hs;
clear WWM NWM JWM L1 sum5 Hs Hm


%---- MAKE TSM OF GRIDPOINT MONTHLY DATA, IN CLIMATE UNITS

%-- Loop over gridpoints
for n = 1:mG;
   i1=1+(n-1)*12;  % starting col of data in GPT1a, and target col in GPT1
   i2=i1+11;  % ending column
   Gtemp=GPT1a(:,i1:i2);  % standardized anomalies, years x months
   Gtempm = repmat(GPTm(n,:),nref,1); % dupe the row of long-term monthly means to a tsm
   Gtemps = repmat(GPTs(n,:),nref,1); % likewise for the standard deviations
   GPT1(:,i1:i2)=Gtempm + (Gtemp .* Gtemps);
end
clear i1 i2 Gtemp Gtempm Gtemps ;
   

%-- ALL-STATIONS ANALYSIS

%---- Associate a A "NEAREST" gridpoint with each station
[Ytemp,inear]=min(Dist1); % inear is rv telling nearest gridpoint to each station
clear Ytempl

%---- Compute unadjusted mean and adjusted mean for each station
%
% strategy
% Pull the monthly data for reference period for each station.
% Make logical pointer matrix to non-nan data in each month
% Store valid sample size for ref period
% Compute the reference-period monthly means for the station
%		also the std devs; store
% Use pointer to loop over months and compute sub-ref-pd monthly
%		means for the master station; same for st devs
% Compute adjustment ratios for means and std devs
% Compute adjusted station means and std devs and store

% Allocate
umean=repmat(NaN,num2,12); % unadjusted long-term monthly means, each station
ustdev=repmat(NaN,num2,12); % unadjusted standard devs
usize=repmat(NaN,num2,12); % sample size (no of non-NaN years for each month of ref pd
amean=repmat(NaN,num2,12); % adjusted long-term means (by ratio method)
astdev=repmat(NaN,num2,12); % adjusted stand devs (by lsq fit of stdev vs adj mean)
ratmn=repmat(NaN,num2,12); % adjustment ratio (ratio of full ref pd to sub-period mean
%    for nearest gridpoint-master series

% Make row pointer to reference-period rows in all-station tsm
igo=pdref(1)-pdstore(1)+1;
isp=pdref(2)-pdstore(1)+1;
irows=(igo:isp);


%-- ALL-STATIONS ANALYSIS: COMPUTE LONG-TERM UNADJUSTED AND ADJUSTED STATISTICS

for nn=1:num2; % loop over all stations 
   jthis=nn; % index to station number
   jcols=J1(jthis,1):J1(jthis,2); % pointers to this station's data in D
   X=	D(irows,jcols); % reference period monthly data
   
   % Unadjusted mean and std dev
   mn1=nanmean(X); % unadjusted ref pd means
   std1=nanstd(X); % unadjusted std devs
   
   % Pointer to NaNs and not-NaNs in station tsm
   L1=~isnan(X); % to not-NaN
   L2=~L1; %  to NaNs
   sum2=sum(sum(L2)); % total number of monthly values NaN 
   sum1=sum(sum(L1)); % total number of not-NaN monthly values
   
   % Store unadjusted means and std devs, and sample size
   umean(nn,:)=mn1;
   ustdev(nn,:)=std1;
   usize(nn,:)=sum(L1);  % store sample size (yrs in ref pd with valid data)
   
   % Identify nearest gridpoint
   iclose = inear(nn);
   
   % Build pointer to cols in GPT1 holding nearest gridpoints monthly data
   jgo1=1+(iclose-1)*12; % starting column
   jsp1=jgo1+11; % ending column
   
   % Pull the gridpoints monthly data
   YYM = GPT1(:,jgo1:jsp1);
   
   % Make NaNs match those months/years in the station series
   % Pull nearest gridpoint's master-interpolated tsm 
   	YYM(L2)=NaN; 
   
   % Compute sub-period means and std devs for gridpoint series
   mny=nanmean(YYM);
	stdy=nanstd(YYM);
   
   % Compute ratio of full-ref-period to sub-period mean at gridpoint 
   % This is the adjustment ratio
   ratmn(nn,:)=GPTm(iclose,:) ./ mny;
   
   % Compute adjusted mean
   amean(nn,:)=ratmn(nn,:) .* mn1;
   
   %-- Use the regression eqns of std vs adjusted mean to get adjusted standard dev
   astdev(nn,:) = sum(B .*  [ones(1,12);  amean(nn,:)]);
   % Above, B is 2 x 12, and amean (nn,:) is 1 x 12     
   % Store unadjusted and adjusted stats, and the adjustment ratios
   UMN(jthis,:)=umean(nn,:);
   USD(jthis,:)=ustdev(nn,:);
   AMN(jthis,:)=amean(nn,:);
   ASD(jthis,:)=astdev(nn,:);
   RATMN(jthis,:)=ratmn(nn,:);
   USIZE(jthis,:)=usize(nn,:);
end


%-- ALL-STATIONS ANALYSIS: COMPUTE TSM OF STANDARDIZED ANOMALIES

% Compute standardized departure time series for stations.
% Fill matrices E and F with standardized departure series. 

for nn=1:num2; % Loop over stations 
   ithis=nn; % rename sequential station number for convenience;
   irows=I(ithis,1):I(ithis,2); % points to rows in D holding this station's data
   jcols=J1(ithis,1):J1(ithis,2); % ...cols in D ....
   X=D(irows,jcols); % Monthly data, entire data record for this station
   [mX,nX]=size(X);
   
   % Uadjusted stats for this station
   mnu =UMN(ithis,:); % mean
   sdu =USD(ithis,:); % std dev
   
   % Adjusted stats for this station 
   mna = AMN(ithis,:); % mean
   sda = ASD(ithis,:); % std dev
   
   % Dupe rvs of means and standard deviations to matrices same size as X
   if kopts(1)==1;  % use adjusted means and std devs
      MNS=repmat(mna,mX,1); % dupe rv to matrix
      SDS=repmat(sda,mX,1);
   elseif kopts(1)==2; % use unadjusted means and std devs
      MNS=repmat(mnu,mX,1);
      SDS=repmat(sdu,mX,1);
   elseif kopts(1)==3; % use adjusted mean and unadjusted standard dev
      MNS=repmat(mna,mX,1);
      SDS=repmat(sdu,mX,1);
   end
   
   % Compute standardized anomalies
   F(irows,jcols)= (X - MNS) ./ SDS;
   
end


%***************  COMPUTE GRIDPOINT-WEIGHTED STDZD ANOMALIES

disp('WEIGHTING STDZD ANOMALIES TO GRIDPOINTS');

% Allocate
nP=mG*12; % # cols in huge grdpnt matx
P=repmat(NaN,nrows,nP); % to hold gridpoint data in clim units (e.g., inches)
Q=repmat(NaN,nrows,nP); % to hold gridpoint data in stdzd units

% Build pointers for each gridpoint's series in P and Q

% What cols in P and Q will the gridpoint monthly data go into?
j1=(1:12:(nP-11)); % start col for point 1, 2,etc
j2=j1+11; % end col
JP1=[j1' j2'];
% First row of JP1 will give first and last storage row for 
% point 1, second row for point 2 etc.  Points numbered as in 
% input coordinate mtx CG

% What cols in P and Q will hold all Jan data, Feb data, etc?
j1=(1:12);
j1=repmat(j1,mG,1);
j2=(1:mG)';
j2=12*(j2-1);
j2=repmat(j2,1,12);
JP2=j1+j2;
% First col of JP2 will point to cols of P with Jan data for all
% mG gridpoints.  Second col to cols of P with Feb data, etc


% Loop over the 12 months of the year
for kmon=1:12;
   
   disp(['   Month ' int2str(kmon)]);
   jp2=JP2(:,kmon); %cv of cols indices to P,Q for this month
   
   % Pull stdzd anoms for this month of yer
   F1=F(:,J2(:,kmon));
	
	% Pull subset of rows of F1 for which not all values are NaN
	L1=~isnan(F1);
	L1a=(any(L1'))';
	F2=F1(L1a,:);

	% Make logical pointer to F2 with data
	L1=~isnan(F2);

	% Compute matrix of weights for this month of year
	kk=[kopts(3) kopts(4)];
	[WW,NW,JW]=wgtdist1(Dist1,L1,dcrit,dsrch,nmax,kk);

	% Loop over the gridpoints
	for kgrd=1:mG;
		if kmode==1;
			disp(['GRIDPOINT # ' int2str(kgrd)]);
		end

		jpthis=jp2(kgrd); % col in P, Q for this month/gridpoint		

		colsWW=JW(kgrd,1):JW(kgrd,2);  % columns of W with this points wts
		W1=WW(:,colsWW);
		
		% Weight the stations
		G1=F2 .* W1;
      
		% Sum over cols to get cv of weighted stdzd anoms this point/month
		g1 = (sumnan(G1'))';

		% Store weighted anomalies
		Q(L1a,jpthis)=g1;
	end
end


%************** GRIDPOINT-FROM-ALL-STATIONS WEIGHTED LONG-TERM MEANS AND STD DEVS
%
disp('COMPUTING GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS FROM ALL STATIONS');

% Allocate
PM=repmat(NaN,mG,12); % to hold means
PS=repmat(NaN,mG,12); % to hold std devs
W5=repmat(NaN,mG,nuse);

% Grab the desired version of means and standard devs
if kopts(1)==1;
   TMN=AMN;
   TSD=ASD;
elseif kopts(1)==2;
   TMN=UMN;
   TSD=USD;
elseif kopts(1)==3;
   TMN=AMN;
   TSD=USD;
else;
   error('invalid kopts(1)');
end

% Mean
L1=~isnan(TMN');
[WWM,NWM,JWM]=wgtdist1(Dist1,L1,dcrit,dsrch,nmax,kk);
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




%************* WEIGHT GRIDPOINT DATA INTO REGIONAL DATA
disp('CALLING GRD2REG.M TO WEIGHT GRIDPOINT DATA TO REGIONAL DATA');
[P,RA,R,RM,RS,JP3,wt2]=grd2reg(IR2,CG,Q,PM,PS,JP2);


% Optionally substitute zero for negative monthly regional values
if kopts(2)==1; % substitute -- most appropriate for ppt data
	L3=R<0;
	numneg=sum(sum(L3));
	if numneg>0;
	 	R(L3)=0;
	end
	L3=P<=0;
	numneg=sum(sum(L3));
	if numneg>0;
		P(L3)=0;
	end
elseif kopts(2)==2; % do not substitute -- most approp. for tmp data
	% no action
else
	error('Invalid setting for kopts(2)');
end

%********************************************************************
disp('CALLING regsub1 TO GET STATION COUNT MATRICES NSTNS AND NS')
nreg=1; % one regions
IR1=(1:nuse)'; % sequence numbers (relative to reduced set of nuse) of stations
[NSTNS,NS,yr3]= regsub1(nuse,nreg,F,J2,IR1,yr2);

disp('here')




%******************** *****************************
disp('BUILDING OUTPUT FILES')

% RECALL THAT

% names1 is row-cell of names of stations used in the analyis (see iuse)
% IR2 tells which gridpoints are included in computing the regional average
% IM is pointer to rows of names1 telling which are master stations
% C2 long lats of stations
% CG long lats of gridpoints
% NS number of good stations matrix, by region
%
% P tsm of monthly gridpoint climate data
% Q tsm ... as standardized anomalies
% PM long term monthly means of gridpoint data
% PS long term standard devs
% RA regional tsm of stdzd anomalies of monthly data
% R  regional tsm monthly data in climate units
% RM regional long-term means -- used to convert RA to R
% RS regional long-term std devs -- used to convert RA to R
% TMN long-term (ref period) monthly means of station data
% TSD long-term   std devs
% W5 matrix of weights used to compute long-term gridpoint
%	  means from long term station means
% wt2 weights of gridpoints to regions
% yr2 year vector for huge matrix
% yr3 year vector for NS, NSTNS after trimming off all-nan rows


% Station list, with sequence number, long, lat, region station in
% for adjustment of mean and standard dev, and number indicating
% if station a master -- and if so, for what region
file3=['stnlst.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
blnk='  ';
for n=1:nuse;
	str1=sprintf('%2.0f',n); % station number
   norig=iuse(n); % station number in original full sequence
   str4=sprintf('(%2.0f)',norig); 
	str2=sprintf('%s',names1{n}); % station code (file prefix)
	str3=sprintf('%6.2f %5.2f ',C2(n,[1 2])); % long, lat
	% master or not
	j=find(IM==n);
	if isempty(j); % not a master
		str5=sprintf('%s',blnk);
	else
		str5=sprintf('%s','*');
		%str5=sprintf('%2.0f',j); % used to want region in parens
	end
	strall=[str1 str4 str5 str3];
	fprintf(fid3,'%s\n',strall);
end
fclose(fid3);


% Gridpoint list:
% gridpoint #;  weight of gridpoint in region
file3=['grdlst.txt'];
pf3=[path5 file3];
fid3=fopen(pf3,'w');
blnk='  ';
for n=1:mG;
	str1=sprintf('%3.0f  ',n); % gridpoint number
	str2=sprintf('%6.2f %5.2f ',CG(n,[1 2])); % long, lat
	% get weight
	w6=wt2(n);
	str4=sprintf('  %6.4f',w6); % weight
	strall=[str1 str2 str4];
	fprintf(fid3,'%s\n',strall);
end
fclose(fid3);

% Listing of stations and weights used to get long term
% gridpoint means and std devs
file3=['wgts.txt'];
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
      thisname=names1{j(k)};
      nctemp=length(thisname);
      bb1=blanks(6);
      if nctemp>6;
         bb1=thisname(1:6);
      else;
         bb1(1:nctemp)=thisname;
      end
      str4=sprintf('  %s',bb1);
		str5=sprintf('(%7.4f)',w6(k));
		strb=[str3 str4 str5];
		fprintf(fid3,'%s\n',strb);
	end
end
fclose(fid3);

%-- .mat storage of "master" series information 
vlist= 'contents of master.mat';
vlist=char(vlist,'% names1 names of stations used in the analysis; master stations are subset of this');
vlist=char(vlist,'% IM row-index to master stations in names1 as stored here');
vlist=char(vlist,'% CG lons,lats of gridpoints');
vlist=char(vlist,'% GPTm gridpoint long-term means interpolated from master stations');
vlist=char(vlist,'% GPTs gridpoint sdevs interpolated from master stations');
vlist=char(vlist,'% GPT1a tsm of standardized anomalies at gridpoints');
vlist=char(vlist,'% GPT1 tsm of master-interpolated gridpoint data in original units');
setmaster=' names1 IM CG GPTm GPTs GPT1a GPT1 vlist';
eval(['save master ' setmaster]);
   
%-- .mat storage of "station" series information 
vlist= 'contents of station.mat';
vlist=char(vlist,'% names1 names of stations used in the analysis');
vlist=char(vlist,'% iuse row-index crossrefering stns in names1 to names1 as input to regcli4.m');
vlist=char(vlist,'% C2 long,lat of stations in names1');
vlist=char(vlist,'% UMN AMN unadjusted means and adjusted means');
vlist=char(vlist,'% USD ASD unadjusted and adjusted standard deviations');
setstn=' names1 IM iuse C2 UMN AMN USD ASD vlist ';
eval(['save master ' setmaster]);

%-- .mat storage of "regional" info
vlist='contents of regional.mat';
vlist=char(vlist,'% CG lons, lats of gridpoints');
vlist=char(vlist,'% IR2 index to rows of CG telling which gridpoints averaged into region');
vlist=char(vlist,'% RA monthly regional standardized anomalies (year vector is yr2)');
vlist=char(vlist,'% yr2 year vector for RA');   
setreg=' CG IR2 RA yr2 ';
eval(['save regional ' setreg]);
   
   
%  .mat file of monthly regional climate series
flout='reg';
% get cols of data for region
jcols=JP3(1,1):JP3(1,2);
Z=R(:,jcols);
X=[yr2 Z];
% output
eval(['save ' path5 flout ' X']);


%-- Build and save gridpoint .mat files interpolated from full set of stations
for n = 1:mG; % loop over gridpoints
	% build file name
	flout=['grd'  int2str(n)];
	% get cols of data for gridpoint
	jcols=JP1(n,1):JP1(n,2);
	Z=P(:,jcols);
	X=[yr2 Z];
	% output
	eval(['save ' path5 flout ' X']);
end


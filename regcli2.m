function regcli2
%
% Regional PPT or TMP series from station data by standardized normal method.
% Produces multiple regional series. One application is to get monthly
% station climate series weighted to tree sites for later use in 
% spatrec1.m.  Another application is to weight station climate data
% to long-lat gridpoints.  Another application is to get regional
% climate series averaged over specified gridpoints.
%
% Weighting implemented so far is  inverse distance or inverse dist squared
% 
% Meko 2-9-97. Revised 1-23-00.
%
%***************************** INPUT FILES *********************************
%
% 1. fls*.txt -- edix-produced file with list of prefixes of filenames
%		of .mat files of the raw monthly climate data.  First line
%		should give path to input .mat files of monthly climate. Line
%		2 should give path to output files created by regcli2.m. For example:
%
%			d:\sprb\ppt\rawmon\
%			d:\outfls\
%			tucs 1
%			bisb 2
%			doug 3
%			...etc
%		The number of filenames, in this file is ns1. The
% 		number to the right of the file name is ignored by the program
% 2. misc*.mat -- miscellaneous program control
%		CG [mG x 2]r long-lat file in mapping units for the mG gridpoints.
%		C2 [ns1 x 2]r long-lat file in mapping units for the ns1 stations
%			Typically produced as ascii out of the
%			database, then loaded and stored in matlab
%		IR1 [? x nreg]i  index to rows of C2 telling which stations
%			are in each region.  Cols padded with zeros to complete mtx
%			IR1 has vital info for the reference period adjustment, because
%			the makeup of master series differs from region to region.
%			Must associate each station with a master series, and so with
%			a region.  Note however, that all ns1 stations are available
%			for the gridpoint weighting, not just the stations in the 
%			gridpoints "region". 
%		IR2 [? x nreg]i index to rows of CG telling which gridpoints are
%			in each region
%		IM (? x nreg)i  index to rows of C2 telling which stations are to
%			be used as master stations for each region.  Cols padded with
%			zeros to complete matrix
%		pdref (1 x 2)i reference period (start, end years) for adjusting
%			means and standard deviations. Assumed the same for all regions
%		pdstore (1 x 2)i  start, end years for storage matrix that must
%			cover at least the period with any valid data at any station.
%			Assumed the same for all regions.
%		kmode (1 x 1)i mode of run.  kmode==1 means diagnostic. kmode==2
%			means skip diagnostic output. This option not operational now.
%			Just set ==2.
%		kopts (1 x 5)i options
%			kopts(1):==1 --> adjust only the means for anomalous sup-period
%						of reference period; do not adjust the std devs
%				    ==2 --> adjust means and std deviations
%			kopts(2):==1 --> subst. zero for negative monthly regional values
%				    ==2 --> accept negative monthly regional values
%					(note: typically will set ==1 for ppt, ==2 for temper.)
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
%						==2 weighting method inverse distance squared
%			kopts(5) ==1 use adjusted means and/or sdevs-- also see kopts(1)
%						==2 use unadjusted
%
%			dcrit -- critical threshold distance to avoid over-weighting site in
%				distance weighting (km)
%			dsrch -- search distance for weighting (km)
%			nmax --  max number of stations to weight to a gridpoint
%
%******************** USER PROMPTED FOR *****************
%			kpass    ==1 "pass 1" mode. Build master series for each region.
%							  In this mode, 
%						==2 "pass 2 mode". Read in master series generated in
%							  in pass 1
%							 files, and read master files off the indicated
%							 .mat files
%			filecode  one-letter code used in output files prefixes to
%					help in documenting run.  For example, might use
%					"m" for a pass 1 run , then "a" for a run using 
%					no mean adjustment and "b" for a run with mean adjustment
%
%			mstpref  string file prefix needed for pass 2 run. The prefix
%					gives first part of .mat file names holding individual
%					regional master series built in earlier pass 1 run. For
%					example, "mstm" would mean the regional files to be
%					read would be mstm1.mat, mstm2.mat ...
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
%******************** DEFINTIONS OF COMPUTED VARIABLES *******************
%
% A distance mtx, gridpoints to stations (km)
% AMN (ns1 x 12) adjusted long-term means, all stations
% ASD (ns1 x 12) adjusted long-term std devs
%- B1 (ns1 x 8)s .mat file prefixes
% mG -- number of gridpoints (or tree locations)
% D huge monthly mtx, original data
% E huge ..., unadjusted standardized anomalies
% F huge ..., 
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

t0=clock;

%
% Allocate
blnks=blanks(8);
B1=blnks(ones(1500,1),:); % Initialize to allow for 200 filenames
sblank=sprintf('%\n'); % used for blank line to screen



% Prompt for whether pass1 or pass 2.
kpass=input('Pass 1 or Pass 2 run (1 or 2)? ');
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
Ltemp1=[exist('IM') exist('IR1') exist('IR2') exist('CG') exist('C2')];
Ltemp2=[exist('pdref') exist('pdstore')];
Ltemp3=[exist('kmode') exist('kopts')];
Ltemp=[Ltemp1 Ltemp2 Ltemp3];
Ltemp=all(Ltemp==1);
if ~Ltemp
	clc
	disp('The following variables should be in the misc.mat file:');
	disp('IM,IR1,IR2,CG,C2,pdref,pdstore,kmode,kopts');
	error('All the above were not in the file');
end
clear Ltemp1 Ltemp2 Ltemp3 Ltemp file2

% Check input arguments


[mG,nG]=size(CG); % x,y coordinates of grids
if nG~=2 ;
	error('CG must have col size 2')
end
if  any(any (CG(:,1)>0));
	error('CG long should be negative')
end
if  any(any(CG(:,2)>90));
	error('CG(:,2) should be all 90 deg or less (latitude)');
end
clear nG

%[mC2,nC2]=size(C2);
%if mC2~=ns1;
%	error('Row size of C2 not equal to number of .mat filenames');
%end
%if nC2~=2;
%	error('C2 should be 2-column');
%end
%if  any(any (C2(:,1)>0));
%	error('C2 long should be negative')
%end
%if  any(any(C2(:,2)>90));
%	error('C2(:,2) should be all 90 deg or less (latitude)');
%end
%clear mC2,nC2


% IM says which stations are to be used in interpolating the
% regional master series 
[mIM,nIM]=size(IM);
nreg=nIM; % number of regions
if nIM~=nreg;
	error('Col size if IM must equal number of regions')
end
if any(any(isnan(IM)));
	error('IM should be zero-padded, not Nan-padded')
end
if any(any(IM>ns1));
	error('An element of IM larger than highest station seq no.');
end
% Compute number of stations in master network
L1=IM>0;
nmast=sum(sum(L1)); % scalar 
if nmast>ns1
	error('nmast greater than number of stations in dataset');
end
clear L1



% On pass 2 IR1 says which stations are considered to be in each
% region for purposes of evaluating station density change over time
% via NS.
%
% On pass 1, IR1 is ignored, and set equal to IM.  This means only
% the master stations are used in interpolating the regional series
% in pass 1.
[mIR1,nIR1]=size(IR1);
if any(any(isnan(IR1)));
	error('IR1 should be zero-padded, not Nan-padded')
end
if any(any(IR1>ns1));
	error('An element of IR1 larger than highest station seq no.');
end

% Number of stations in each region
L1=IR1~=0;
ns2=sum(L1);
if kpass==1; % pass 1
	IR1=IM;
	[mIR1,nIR1]=size(IM);
	L1=IR1>0;
	ns2=sum(L1);
end

% IR2 says which gridpoints are in which region. Info needed for
% weighting gridpoint data in to regional series
[mIR2,nIR2]=size(IR2);
if any(any(isnan(IR2)));
	error('IR2 should be zero-padded, not Nan-padded')
end
if any(any(IR2>mG));
	error('An element of IR2 larger than highest grdpnt seq no.');
end
if nIR2~=nreg,
	error('IR2 should be same col size as IR1')
end
% Compute number of gridpoints in each region
L1=IR2>0;
npts1=sum(L1); % rv 
if npts1>mG
	error('npts1 greater than mG -- row size of CG');
end
clear L1





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
if kopts(1)<1 | kopts(1)>2 
	error('Illegal setting of kopts(1)')
end
if kopts(2)<1 | kopts(2)>2;
	error('Illegal setting of kopts(2)')
end
if kopts(3)<1 | kopts(3)>2;
	error('Illegal setting of kopts(3)')
end
if kopts(4)~=1 & kopts(4)~=2;
	error('Illegal setting of kopts(4)')
end
if kopts(5)<1 | kopts(5)>2;
	error('Illegal setting of kopts(5)')
end

if kpass==1; % "pass 1" mode. Will build master series.
	regpref='mst'; % prefix to output "regional" .mat files
	% Re-set the option for whether to use unadjusted or adjusted
	% means and/or standard deviations.  In pass 1 mode, want to
	% use unadjusted means and standard devs always
	kopts(5)=2; 

elseif kpass==2; %"pass 2" mode. Already built master in a 
	% pass 1 run. Now building regional series using full data.
	% read prefix of master series (e.g., mstm -- which would
	% mean master series are in mstm1.mat, mstm2.mat, ...
	mstpref=input('Master series prefix (e.g., mstm): ','s');
	regpref='reg';
else
	error('kpass must be 1 or 2');
end

if kpass==1;
	disp('YOU ASKED FOR A PASS 1 RUN');
else;
	disp('YOU ASKED FOR A PASS 2 RUN');
end

if kopts(5)==2;
	disp('COMPUTATIONS WILL USE UNADJUSTED MEANS AND STD DEVS');
else
	if kopts(1)==1;
		disp('COMPUTATIONS USE ADJUSTED MEANS, BUT UNADJUSTED STD DEVS');
	else
		disp('COMPUTATION USE ADJUSTED MEANS AND STD DEVS');
	end
end
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
	
	disp('Number of stations in each region:')
	disp(ns2)
	
	disp('Number of gridpoints in each region: ')
	disp(npts1)
	
	disp('Number of stations in in master network:')
	disp(nmast);
	
	disp('Reference period for master series')
	disp([pdref(1) pdref(2)]);
	
	disp('Period for huge matrix');
	disp([pdstore(1) pdstore(2)]);
	
	disp('Settings for call to wgtdista.m')
	disp([' dcrit = ' num2str(dcrit) ' km'])
	disp([' dsrch= ' num2str(dsrch) ' km'])
	disp([' nmax= ' int2str(nmax) ' stations'])
	disp([' Option for ignoring search radius: ' int2str(kopts(3))])
	disp([' Option for weighting method: ' int2str(kopts(4))])
	disp([' Option for using adjusted vs unadjusted: ' int2str(kopts(5))])
	
	disp('Press any key to continue')
   pause
   clc
end

% Allocate
a=NaN;
AMN=a(ones(ns1,1),ones(12,1));
ASD=a(ones(ns1,1),ones(12,1));
UMN=a(ones(ns1,1),ones(12,1));
USD=a(ones(ns1,1),ones(12,1));
RATMN=a(ones(ns1,1),ones(12,1));
RATSD=a(ones(ns1,1),ones(12,1));
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
F=a(ones(nrows,1),ones(ncols,1));  % monthly data, adjusted anoms
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
	clear Z;
	fln=strtok(B1(n,:));
   eval(['DATA=load(''' path3 fln ''');']);
   LZ=isfield(DATA,'Z');
   LX=isfield(DATA,'X');
   LY=isfield(DATA,'Y');
   Lxyz=[LX LY LZ];
   if sum(Lxyz)==1;
   elseif sum(Lxyz)==0;
      error([fln ' does not contain X Y or Z']);
   else;
      error([fln ' contains more than one of variables X , Y, Z']);
   end;
   if LZ;
      Z=DATA.Z;
   elseif LY;
      Z=DATA.Y;
   elseif LX;
      Z=DATA.X;
   else;
      error('No DATA.X, Y or Z');
   end;
   clear DATA;
   
   
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
if kpass==1;
	disp('A PASS 1 RUN -- OBJECTIVE IS TO BUILD MASTER SERIES');
	disp('COMPUTING STDZD ANOMALY SERIES FOR STATIONS');
else
	disp('A PASS 2 RUN -- OBJECTIVE IS GRIDPOINT AND REGIONAL SERIES');
	disp('READING MASTER SERIES & COMPUTING STDZD ANOMALY SERIES');
end

for n=1:nreg; % loop over regions
	disp(['  Region # ' int2str(n)]);
	% get index to stations for master
	im=IM(:,n);
	im(im==0)=[];
	% get index to stations in region
	ir = IR1(:,n);
	ir(ir==0)=[];

	% Calc number of stations in regional series and in master
	num2=length(ir); % # stns in region
	num3=length(im); % # stns in master
	B2=B1(im,:); % stations for master
	B3=B1(ir,:); % stations in region 
	if kmode==1;
		
		disp('Series specified for master')
		disp(B2);
		
		disp('Series specified for region')
		disp(B3);
		
		disp('Press any key to continue')
		pause	
	end

	% size and allocate
	numcol=num3*12; % # cols in master storage mtx
	M=a(ones(nref,1),ones(numcol,1)); %  ref pd raw data, years by months
	LL=a(ones(nref,1),ones(numcol,1)); %  ref pd missing value pointer
	V=a(ones(nref,1),ones(numcol,1)); % ref period station stdzd anomalies
	VM=a(ones(nref,1),ones(12,1)); % mean stdzd anomaly series for the master
	YM=a(ones(nref,1),ones(12,1)); % VM converted to original climate units
	Mrefm=a(ones(num3,1),ones(12,1)); % ref period long-term monthly means
	Mrefs=a(ones(num3,1),ones(12,1)); % ref period long term std devs
	gm=a(:,ones(12,1)); % global mean monthly ppt (average over stations)
	gs=a(:,ones(12,1)); % global mean std dev


	if kpass==1; % pass 1 mode
		% Get the monthly data for the master stations and store in matrix M
		Lmast = yr2>=pdref(1) & yr2<=pdref(2); % ref pd rows in the mtx
		for nn = 1:num3;
			j1=J1(im(nn),:); % start,end cols for this station in the huge mtx
			kcols=j1(1):j1(2);
			Z=D(Lmast,kcols); % the monthly climate data for ref period
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
		end

		%  Average stdzd anomalies over master stations
		for k = 1:12; % loop over the 12 months of the year
			jj1=(1:12);
			JJ1=jj1(ones(num3,1),:); 
			jj2=(1:num3)';
			jj3=12*(jj2-1);
			JJ3=jj3(:,ones(12,1));
			JJ4=JJ1+JJ3;
			jcols=(JJ4(:,k))';
			H=V(:,jcols);
			vm=(meannan(H'))'; % vector of mean stdzd anoms, this month of yr
			VM(:,k)=vm; % slap this month's onto the matrix
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
	else; % kopt(6)==2, and want to read in master series and
		% compute ref period means and standard deviations
		clear Z;
		flmst=[path5 mstpref int2str(n)];
		eval(['load ' flmst]); % data in Z
		if (exist('Z')~=1)
			disp(['Region ' int2str(n)]);
			error(['No Z in: ' flmst]);
		end

		% Make year vector; store first year and last year; calc and store
		% row reference info
		yr5 = Z(:,1);
		nyrs5=length(yr5);
		if yr5(1)>pdref(1) | yr5(nyrs5)<pdref(2)
			disp(flmst)
			error('Master does not cover reference period');
		end
		d1=diff(yr5);
		if ~all(d1==1);
			disp(flmst);
			error('Some internal year missing from master');
		end

		% Grab reference period years of master
		Ltemp=yr5>=pdref(1) & yr5<=pdref(2);
		YM=Z(Ltemp,2:13);

		% Check that no NaNs in reference master
		Ltemp=isnan(YM);
		if any(any(Ltemp));
			disp(flmst);
			error('NaN in the master in reference period');
		end

		% Compute monthly reference period means and std devs for master
		gm = mean(YM);
		gsnew=std(YM);
	end

	%  MASTER SERIES AND ASSOCIATED VARIABLES FOR THIS REGION BUILT
	% YM tsm of regional climate variable for the master series, in original
	% 		units,for the reference period
	% gm rv of master-series monthly means for reference period
	% gsnew rv of master-series monthly std devs for ref pd


	%
	%	******************************************************************
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

	umean=a(ones(num2,1),ones(12,1));
	ustdev=a(ones(num2,1),ones(12,1));
	usize=a(ones(num2,1),ones(12,1));
	amean=a(ones(num2,1),ones(12,1));
	astdev=a(ones(num2,1),ones(12,1));
	ratmn=a(ones(num2,1),ones(12,1));
	ratsd=a(ones(num2,1),ones(12,1));

	igo=pdref(1)-pdstore(1)+1;
	isp=pdref(2)-pdstore(1)+1;
	irows=(igo:isp);

	for nn=1:num2; % loop over all stations in the region
		jthis=ir(nn); % index to station number
		jcols=J1(jthis,1):J1(jthis,2);
		X=	D(irows,jcols); % reference period monthly data
		L1=~isnan(X);
		usize(nn,:)=sum(L1);  % store sample size (yrs in ref pd with valid data)
		[mn1,msize]=meannan(X); % unadjusted ref pd means; second arg not used
		[std1,msize]=stdnan(X); % unadjusted std devs
		umean(nn,:)=mn1;
		ustdev(nn,:)=std1;


		% Master series means and std devs for the sub period of ref period
		L2=~L1; % pointer to NaNs in station series
		sum1=sum(sum(L2)); 
		YYM=YM;
		YYM(L2)=a(ones(sum1,1),:); 
		[mny,msize]=meannan(YYM);
		[stdy,msize]=stdnan(YYM);

		% Ratios for adjusting station means and std devs
		ratmn(nn,:)=gm ./ mny;
		ratsd(nn,:)=gsnew ./ stdy;

		% Compute adjusted means and std devs
		amean(nn,:)=ratmn(nn,:) .* mn1;
		astdev(nn,:)=ratsd(nn,:) .* std1;

		% Store unadjusted and adjusted stats, and the adjustment ratios
		UMN(jthis,:)=umean(nn,:);
		USD(jthis,:)=ustdev(nn,:);
		AMN(jthis,:)=amean(nn,:);
		ASD(jthis,:)=astdev(nn,:);
		RATMN(jthis,:)=ratmn(nn,:);
		RATSD(jthis,:)=ratsd(nn,:);
		USIZE(jthis,:)=usize(nn,:);
	end


	%********************************************************************
	%
	% Compute and store means and standard deviations averaged over
	% all stations marked by IR1 in region 1, region 2, etc. These
	% arithmetic means are not used elsewhere in the program, but 
	% might be interesting to compare with the regional monthly means
	% and standard devs by the distance weighting and latitude weighting
	% algorithm
	ustats=a(ones(2,1),ones(12,1)); % storage for means in row 1, sdevs in 2
	astats=a(ones(2,1),ones(12,1)); % storage for means in row 1, sdevs in 2
	% Unadjusted 
	ustats(1,:)=mean(umean);
	ustats(2,:)=mean(ustdev);
	astats(1,:)=mean(amean);
	astats(2,:)=mean(astdev);



	%*******************************************************************
	%
	% Compute standardized departure time series for stations.
	% Fill matrices E and F with standardized departure series. 

	for nn=1:num2; % Loop over stations in region
		ithis=ir(nn);
		irows=I(ithis,1):I(ithis,2); % points to rows in D holding this station's data
		jcols=J1(ithis,1):J1(ithis,2); % ...cols in D ....
		X=D(irows,jcols); % Monthly data, entire data record for this station
		[mX,nX]=size(X);

		% Using unadjusted stats
		mns = umean(nn,:);
		sds = ustdev(nn,:);
		MNS=mns(ones(mX,1),:); % rv to matrix
		SDS=sds(ones(mX,1),:);
		E(irows,jcols)= (X - MNS) ./ SDS;

		% Using adjusted stats
		mns = amean(nn,:);
		if kopts(1)==1; % use adjusted mean, but not unadjusted st dev
			sds=ustdev(nn,:);
		else; % use adjusted mean and adjusted std dev
			sds = astdev(nn,:);
		end
		MNS=mns(ones(mX,1),:); % rv to matrix
		SDS=sds(ones(mX,1),:);
		F(irows,jcols)= (X - MNS) ./ SDS;
	end

end; % of do n= over regions


%
%
%***************  COMPUTE GRIDPOINT-WEIGHTED STDZD ANOMALIES

disp('WEIGHTING STDZD ANOMALIES TO GRIDPOINTS');

% Allocate
nP=mG*12; % # cols in huge grdpnt matx
P=a(ones(nrows,1),ones(nP,1)); % to hold gridpoint data in clim units
Q=a(ones(nrows,1),ones(nP,1)); % to hold gridpoint data in stdzd units

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
j1=j1(ones(mG,1),:);
j2=(1:mG)';
j2=12*(j2-1);
j2=j2(:,ones(12,1));
JP2=j1+j2;
% First col of JP2 will point to cols of P with Jan data for all
% mG gridpoints.  Second col to cols of P with Feb data, etc

% Compute distance matrix (km) from gridpoints to stations
A=gcdist(CG,C2);

% Loop over the 12 months of the year
for kmon=1:12;

	disp(['   Month ' int2str(kmon)]);
	jp2=JP2(:,kmon); %cv of cols indices to P,Q for this month

	% Pull the subset of station stdzd anomalies, all yr in E or F
	if kopts(5)==1; % use unadjusted means and std devs
		F1=E(:,J2(:,kmon));
	else; % use adjusted means and/or std devs
		F1=F(:,J2(:,kmon));
	end

	% Pull subset of rows of F1 for which not all values are NaN
	L1=~isnan(F1);
	L1a=(any(L1'))';
	F2=F1(L1a,:);

	% Make logical pointer to F2 with data
	L1=~isnan(F2);

	% Compute matrix of weights for this month of year
   kk=[kopts(3) kopts(4)];
   
     
   
	[WW,NW,JW]=wgtdist1(A,L1,dcrit,dsrch,nmax,kk);

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

%************** GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS
%
disp('COMPUTING GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS');

% Allocate
PM=a(ones(mG,1),ones(12,1)); % to hold means
PS=a(ones(mG,1),ones(12,1)); % to hold std devs
W5=a(ones(mG,1),ones(ns1,1));

% Optionally use adjusted or unadjusted means and/or std devs
if kopts(5)==1;
	TMN=AMN;
	if kopts(1)==1; % use adjusted means, but unadjusted std devs
		TSD=USD;
	else; % use adjusted std dev as well
		TSD=ASD;
	end
else
	TMN=UMN;
	TSD=USD;
end


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


%************* WEIGHT GRIDPOINT DATA INTO REGIONAL DATA
disp('CALLING GRD2REG.M TO WEIGHT GRIDPOINT DATA TO REGIONAL DATA');
[P,RA,R,RM,RS,JP3,wt2]=grd2reg(IR2,CG,Q,PM,PS,JP2);


% Optionally substitute zero for negative monthly regional values
if kopts(2)==1; % substitute -- most appropriate for ppt data
	L3=R<0;
	numneg=sum(sum(L3));
	if numneg>0;
	 	R(L3)=zeros(numneg,1);
	end
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

%********************************************************************
disp('CALLING regsub1 TO GET STATION COUNT MATRICES NSTNS AND NS')
[NSTNS,NS,yr3]= regsub1(ns1,nreg,F,J2,IR1,yr2);


%******************** *****************************
disp('BUILDING OUTPUT FILES')
 
%****** OUTPUT THE RESULTS **************************
% This section written and tested outside of regcli2.out, then
% added later. 
% 
% RECALL THAT

% IR1 tells which stations in which region 
% IR2 tells which gridpoints in which region
% B1 station filename list
% C2 long lats of stations
% CG long lats of gridpoints
% NSTNS number of good stations matrix, all regions combined
% NS number of good stations matrix, by region

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

% filecode -- 1-char letter added to file prefixes for record keeping



% Reduce row size of P, Q, R, RA by lopping off years for which
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
etime(clock,t0)
fclose all;


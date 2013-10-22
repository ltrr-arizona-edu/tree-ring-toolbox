function regcli5
%
% Interpolate monthly PPT or TMP series from station data by standardized
% anomalies method. Adapted from regcli2.m to work with
% developing site-interpolated climate series for ~800 tree ring sites by
% ~1500 GHCN climate stations. One application is to get monthly
% station climate series weighted to tree sites for later use in 
% spatrec1.m.  Another application is to weight station climate data
% to long-lat gridpoints. Adapted from regcli2.m.  Difference is that 
% regcli5.m only interpolates to gridpoints or sites, and does not take
% the additional step of averaging over points or sites to create regional
% series.
%
% Weighting implemented so far is  inverse distance or inverse dist squared
% 
% Meko 4-27-00
%
%***************************** INPUT FILES *********************************
%
% 1. fls*.txt -- edix-produced file with list of prefixes of filenames
%		of .mat files of the raw monthly climate data.  First line
%		should give path to input .mat files of monthly climate. Line
%		2 should give path to output files created by regcli5.m. For example:
%
%			c:\data\ghcn\pcpv2\matfls\, or ...\tmpv2e\matfls\, or ...tmpv2d\matfls\
%			c:\wrk0\outfls1\,      or  c:\wrk0\outflse\  or c:\wrk0\outflsd\
%			P1   or D1, or E1 
%			P3
%			P8
%			...etc
%		The number of filenames, in this file is ns1.  The file and directory
%		names will vary depending on the data type.  For example, the above are for
%		precip.  But for maximum temp, filenames start with D instead of P.
% 2. misc*.mat -- miscellaneous program control *=P,D,E for precip, max T, min T
%		CG [mG x 2]r long-lat file in mapping units for the mG gridpoints.
%		C2 [ns1 x 2]r long-lat file in mapping units for the ns1 stations
%			Typically produced as ascii out of the
%			database, then loaded and stored in matlab
%		pdstore (1 x 2)i  start, end years for internal storage matrix that must
%			cover at least the period with any valid data at any station
%		kmode (1 x 1)i mode of run.  kmode==1 means diagnostic. kmode==2
%			means skip diagnostic output. This option not operational now.
%			Just set ==2.
%		kopts (1 x 5)i options
%			kopts(1):==1 --> Obsolete, set ==1.
%			kopts(2):==1 --> subst. zero for negative monthly regional values
%				      ==2 --> accept negative monthly regional values
%					(note: typically will set ==1 for ppt, ==2 for temper.)
%			kopts(3) ==1 nearest nmax stns, disregarding search radius is to
%							be used in weighting gridpoints.  Will always have
%							nmax stations in weighting, unless fewer than nmax
%							are available in that month/year for the whole ns1
%							station set.  Every gridpoint will have the same 
%							time coverage in its .mat output file, and this coverage
%							will be from the earliest data at any station to the
%							latest data at any station.
%				    	==2 nearest up-to-nmax stns,  and any station used in 
%							the weighting must be in the search radius. See 
%							"search radius" in NOTES.  
%							all gridpoints can be estimated
%			kopts(4) ==1 weighting method inverse distance
%						==2 weighting method inverse distance squared
%			kopts(5) ==2 Obsolete-- set ==2;
%			dcrit -- critical threshold distance to avoid over-weighting site in
%				distance weighting (km)
%			dsrch -- search distance for weighting (km)
%			nmax --  max number of stations to weight to a gridpoint
%
%******************** USER PROMPTED FOR *****************
%			filecode  one-letter code used in output files prefixes to
%					help label series. Suggest "P" for precip, "D" for max T,
%              'E' for min T, "M" for mean T. Letter will be first char in
%               output gridpoint or site-interpolated .mat filenames
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
%		(* ) if used () if not
% Smn.txt ----- matrix of station means
% Sstd.txt --- stn seq#, matrix of station standard dev
% Nref.txt  -- stn seq#, 12 sample sizes -- number of years ref-pd means
%			and standard deviations based on
%
% GRIDPOINT INFO
%
% Gmn.txt gridpoint interpolated monthly means
% Gsd.txt gridpoint interpolated std devs
% Glist.txt -- point #, lon, lat (5 station seq nos)-- these for means and std devs
% Gwt.txt point#, lon,lat, (5 distances (km)), (5 weights)-- these on means and std devs
%
% Gcover.mat quality of data coverage for each gridpoint
%  Lsome logical vector, length mG.
%		==1  has at one month/year with a predictor station in search radius
%     ==0  does not have any months or years with a station in the search radius
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
% The station means and std devs are computed on all available data.
%
% The standardized anomaly series are then weighted over stations around a 
% gridpoint  (subject to constraint on search radius) to produce a gridpoint
% average standardized anomaly monthly time series.
%
% The gridpoint standardized anomalies are converted to original units
% (e.g., in for ppt) using a weighted average of reference-period means
% and standard deviations for stations in the search radius, dsrch, around
% the gridpoint.
%
% "Gridpoints" might be long-lat points for a regular grid, or might be
% tree locations if want site-weighted climate series. In this case, the 
% gridpoints are irregularly spaced.
%
% Search radius.  If kopts(3)==2, and for some month/year no stations are in
% the search radius, the interpolated value is NaN.  If no stations are in 
% the search radius in any month/year, a .mat file is created with Z==[].
% The output variable Lgood in Gcover.mat flags stations with some data as 1,
% and NaN stations as 0.
%
% Output gridpoint .mat files.  Names G?nnn.mat, where ? indicates data type. 
% Thus GP121.mat is precip for gridpoint 121.   The time coverage by the .mat 
% files is uniformly pdstore if kopts(3)==1.  But if kopts(3)==2, the .mat file
% is trimmed to remove leading and trailing all-NaN years.  Internal NaN
% years or months are retained, however.
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
% ncols -- number of cols in huge monthly multi-station matrix
% npts1 (rv)  -- number of gridpoints in each region
% nref -- # yr in reference period
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
% W5 -- weights on stations to gridpoints for long-term means
% yr2 -- year vector for huge multistation monthly matrix

t0=clock;

% Type of climate variable
kmen1=menu('Choose type of climate variable',...
   'P == precip',...
   'D == maximum Temperature',...
   'E == minimum Temperature',...
   'M == mean Temperature');
if kmen1==1; 
   ctype='P';
   inpmat='c:\wrk0\outfls1\';
   stnfile='stninf.mat';
   sitefile='siteinf.mat';
   lookfile='looky.mat';
elseif kmen1==2;
   ctype='D';
   inpmat='c:\wrk0\outflsd\';
   stnfile='stninfd.mat';
   sitefile='siteinfd.mat';
   lookfile='lookyd.mat';
elseif kmen1==3;
   ctype='E';
   inpmat='c:\wrk0\outflse\';
   stnfile='stninfd.mat'; % because station info same for d,e,m
   sitefile='siteinfe.mat';
   lookfile='lookye.mat';
elseif kmen1==4;
   ctype='M';
   inpmat='c:\wrk0\outflsm\';
   stnfile='stninfm.mat'; % 
   sitefile='siteinfm.mat';
   lookfile='lookym.mat';
end;

% Check that stnfile exists in cwdir -- sequence of files ending with 
% pull39a.m and prepfls1.m should have been run
if ~exist(stnfile,'file')==2;
   error([stnfile ' does not exist in current working directory']);
end;

% Allocate
blnks=blanks(8);
B1=blnks(ones(1500,1),:); % Initialize to allow for 1500 filenames
sblank=sprintf('%\n'); % used for blank line to screen

% Prompt for diagnostics mode
kdiag=questdlg('Diagnostic Mode?');
if strcmp(kdiag,'Yes');
   if exist(lookfile,'file')==2;
      eval(['delete ' lookfile ';']);
   end;
end;

% Build 2-letter prefix start for output .mat files.  
filecode=['G' ctype];   % G for point, ctype for data type


%********************************************************************
% Read list of station .mat files; count # of files; store file names
%First line of the list is the path to input files
%Second line is the path to output files
[file1,path1]=uigetfile(['fls' ctype '.txt'],'List of input .mat files');
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

% Build cell matrix of strings that are prefixes of .mat station monthly 
% climate files
disp('Storing cell matrix of stn filenames');
flnm=fliplr(deblank(fliplr(B1))); % Make sure no leading blanks
flnm=cellstr(flnm); % to cell.  cellstr removes any trailing blanks

if strcmp(kdiag,'Yes');
   vlist='B1 -- string matrix of station filename prefixes';
   vlist=char(vlist,'ns1 -- number of stations, or rows, in B1');
   vlist=char(vlist,'path3 -- to directory with input .mat station data files');
   vlist=char(vlist,'path5 -- to directory with output .mat point or grid .mat files');
   eval(['save ' lookfile ' vlist B1 ns1 path3 path5;']);
else;
   %save stninf flnm -append;
   eval(['save ' stnfile ' flnm  -append;']);
   clear flnm;
end;

% B1 now contains file name prefixes
% ns1 is the number of files
% path3 is the path to the .mat input station data files


% Get miscellaneous program control
[file2,path2]=uigetfile(['misc' ctype '.mat'],'Miscellaneous program control');
pf2=[path2 file2];
eval(['load ' pf2]);


%*********** QC the input args

% Check that all variables wanted from the misc file exist
Ltemp1=[exist('CG') exist('C2')];
Ltemp2=[exist('pdstore')];
Ltemp3=[exist('kmode') exist('kopts')];
Ltemp=[Ltemp1 Ltemp2 Ltemp3];
Ltemp=all(Ltemp==1);
if ~Ltemp
	clc
	disp(['The following variables should be in the misc' ctype '.mat file:']);
	disp('CG,C2,pdstore,kmode,kopts');
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


% Store lon,lat for stations
xy=C2;
if strcmp(kdiag,'Yes');
else;
   %save stninf xy -append;
   eval(['save ' stnfile ' xy -append;']);
   clear xy;
end;

% Prompt for start and end gridpoint for this run
prompt={'Enter starting and ending sequence numbers'};
def={num2str([1 mG])};
dlgTitle='Starting and ending gridpoint or site numbers for this run';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
nboth=str2num(answer{1});
nstart=nboth(1); nend=nboth(2);

if strcmp(kdiag,'Yes') & (nend-nstart)~=0;
   error('Diagnostics mode only valid if one point to be analyzed');
end;
   

if nend>mG;
   error('Specified ending gridpoint number too high');
end;
if nstart>nend;
   error('Starting gridpoint higher than ending');
end;
npnts=nend-nstart+1; % number of points to analyze in this run


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
if kopts(1) ~=1; 
	error('Illegal setting: kopts(1) must set ==1')
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
if kopts(5) ~= 2;
	error('kopts(5) must be set ==2')
end


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
	disp([int2str(mG) ' gridpoints'])
			
	disp('Period for huge storage matrix');
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
UMN=repmat(NaN,ns1,12);
USD=repmat(NaN,ns1,12);
USIZE=repmat(NaN,ns1,12);

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

%*******************************************************************
%
% Compute standardized departure time series for stations.
% Fill matrices E and F with standardized departure series. 

for nn=1:ns1; % Loop over stations in region
   ithis=nn;
   irows=I(ithis,1):I(ithis,2); % points to rows in D holding this station's data
   jcols=J1(ithis,1):J1(ithis,2); % ...cols in D ....
   X=D(irows,jcols); % Monthly data, entire data record for this station
   [mX,nX]=size(X);
   
   % Compute stats
   Lok=~isnan(X);
   USIZE(nn,:)=sum(Lok); % sample size on which mean for each month based on
   mns = nanmean(X);
   UMN(nn,:)=mns; % store with unadjusted means of other stations
   sds = nanstd(X);
   if any(sds==0);
      Lwarn = find(sds);
      
      disp('here');
   end;
   
   USD(nn,:)=sds; % store with unadjusted sdevs ...
   MNS=repmat(mns,mX,1); % rv to matrix
   SDS=repmat(sds,mX,1);
   
   % Handle the case of all years values for a particular month being zero.
   % This gives zero mean and zero standard deviation, and 0/0 ==NaN as
   % the standardized anomaly. Make so that the standardized anomaly is 0.
   Lflag1= X==0 & MNS==0 & SDS==0;  
   Lflag2= X~=0 & MNS~=0 & SDS==0;
   if any(any(Lflag2));
      error('Some month/year has nonzero value but std dev is zero');
   end;
   SDST=SDS;
   if any(any(Lflag1));
      SDST(Lflag1)=1; % set standard deviation to 1 temporarily to avoid div by zero
   end;
   
   E(irows,jcols)= (X - MNS) ./ SDST; % standardized anomalies for indiv stations
end


% BEGIN LOOPING OVER GRIDPOINTS
%***************  COMPUTE GRIDPOINT-WEIGHTED STDZD ANOMALIES

if strcmp(kdiag,'Yes');
   N12 = repmat(NaN,nrows,12); % to hold diagnostics mode sample size 
end;

disp('WEIGHTING STDZD ANOMALIES TO GRIDPOINTS');

for nnn=nstart:nend; % loop over gridpoints or target sites
   disp(['  Gridpoint # ' int2str(nnn)]);
   % Allocate
   nP=1*12; % # cols in grdpnt matx, one col for each month
   P=repmat(NaN,nrows,12); % to hold gridpoint data in original clim units
   	% The nrows rows each are a year.
   Q=repmat(NaN,nrows,12); % to hold gridpoint data in stdzd units
   
   % Build pointers for each gridpoint's series in P and Q
   
   % What cols in P and Q will the gridpoint monthly data go into?
   j1=1; % start col, Jan
   j2=12;   % end col, Dec
   JP1=[j1 j2];
      
   % What cols in P and Q will hold all Jan data, Feb data, etc?
   JP2=1:12;
   % First col of JP2 will point to cols of P with Jan data for the gridpoint
   % Second col to Feb data, etc
   
   % Compute rv of distance matrix (km) from gridpoint to stations
   A=gcdist(CG(nnn,:),C2);
   
   % Set flag if restricting search to only those stations in search radius dsrch, and if
   % have no stations within dsrch of the site or gridpoint
   if all(A>dsrch & kopts(3)==2);
      inone=1;
      if strcmp(kdiag,'Yes');
         error('No stations in search radius. Diagnostic mode not allowed for this point');
      end;
   else;
      inone=0;
   end;
      
   if strcmp(kdiag,'Yes');
      vlist=char(vlist,'A - distance to stations');
      vlist=char(vlist,'nstart - gridpoint number');
      eval(['save ' lookfile '  vlist A nstart -append;']);
   else;
      % Store info on nearest nmax stations and their distances
      if nnn==1; % if first gridpoint
         %if exist('siteinf.mat') ~=2;
         if exist(sitefile) ~=2;
            error([sitefile ' must exist. First build with sitebld.m']);
         end;
         % Initialize matrices to store distance to nearest up-to-nmax stations
         % within search radius dsrch,
         % which stations those are, and how many of those there are
         wlist='site -- structure of site info, from sitebld.m';
         wlist=char(wlist,'Gnm -- cell mtx of prefixes of .mat stge files of monthly data');
         wlist=char(wlist,'Inear -- sequence numbers of nearest nmax stations');
         wlist=char(wlist,'Dnear -- how near (km) are stations indexed in Inear');
         wlist=char(wlist,'Nnear -- how many stations are within dsrch km of point');
         Dnear=repmat(NaN,mG,nmax);
         Inear=repmat(NaN,mG,nmax);
         Nnear=repmat(NaN,mG,1);
         Gnm=cell(mG,1); % to hold prefixes filenames for site-interp series
         eval(['save ' sitefile ' Dnear Inear Nnear Gnm wlist -append;']);
      end;  % if nnn==1
      eval(['load ' sitefile ' Dnear Inear Nnear;']);
      if inone==0; % if any stations within search radius
         [dnear,inear]=sort(A);
         dnear=dnear(1:nmax); inear=inear(1:nmax);
         ikill=dnear>dsrch;
         if any(ikill);
            dnear(ikill)=NaN;
            inear(ikill)=NaN;
         end;
         if all(isnan(dnear));
            nnear=0;
         else;
            nnear=sum(~isnan(dnear));
         end;
      else;
         nnear=0;
         inear=repmat(NaN,1,nmax);
         dnear=repmat(NaN,1,nmax);
      end;
      Nnear(nnn)=nnear;
      Inear(nnn,:)=inear;
      Dnear(nnn,:)=dnear;
      eval(['save ' sitefile ' Dnear Inear Nnear -append;']);
      clear Dnear Inear Nnear dnear inear nnear ikill;
  end;
   
      
   % Loop over the 12 months of the year, if any stations in search radius
   if inone==0;
      for kmon=1:12;
         disp(['   Month ' int2str(kmon)]);
         jp2=JP2(kmon); % of col index to P,Q for this month
         
         % Pull the subset of station stdzd anomalies, all yr in E 
         F1=E(:,J2(:,kmon));
         
         % Pull subset of rows of F1 for which not all values are NaN
         L1=~isnan(F1);
         L1a=(any(L1'))';
         F2=F1(L1a,:);
         if size(F2,1)==1;
            error('Only one year (row of F1) for which not all values (stations)  NaN');
         end;
         
         % Make logical pointer to F2 with data
         L1=~isnan(F2);
         
         % Compute matrix of weights for this month of year
         kk=[kopts(3) kopts(4)]; % options for call to wgtdist1
         [WW,NW,JW]=wgtdist1(A,L1,dcrit,dsrch,nmax,kk);
         
         jpthis=jp2; % col in P, Q for this month/gridpoint		
         
         colsWW=JW(1):JW(2);  % columns of W with this points wts
         W1=WW(:,colsWW);
         
         % Weight the stations
         G1=F2 .* W1;
         
         % Sum over cols to get cv of weighted stdzd anoms this point/month
         g1 = (sumnan(G1'))';
         
         % Change elements of g1 based on zero number of stations from 0 to NaN
         g1(NW==0)=NaN;
         
         % Store weighted anomalies
         Q(L1a,jpthis)=g1;
         
         if strcmp(kdiag,'Yes');
            N12(L1a,kmon)=NW;
            wmonth=repmat(NaN,nrows,ns1);
            wmonth(L1a,:)=WW;
            vlist=char(vlist,'w? -- distance matrix for month ?');
            eval(['w' int2str(kmon) ' = wmonth;']);
            %eval(['save looky  vlist w' int2str(kmon) ' -append;']);
            eval(['save ' lookfile ' vlist w' int2str(kmon) ' -append;']);
            eval(['clear w' int2str(kmon)  ';']);
         end;
         
      end; % loop over months
   else; % No stations are in search radius.  inone is not == 0
      % minor action here, but recall that Q is still all NaNs, with dimensions nrows by 12
      WW=repmat(NaN,nrows,ns1);
      NW=zeros(nrows,1);
      JW=[1 ns1];
      
   end; % if inone==0
   
   
   if strcmp(kdiag,'Yes');
      vlist=char(vlist,'N12 -- sample size (number of predictor statios) in each year');
      vlist=char(vlist,'yr2 -- vector years for N12 and other variables');
      eval(['save ' lookfile ' vlist N12 yr2 -append;']);
      clear N12;
   end;
   
   
   %************** GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS
   %
   disp('COMPUTING GRIDPOINT WEIGHTED LONG-TERM MEANS AND STD DEVS');
   
   % Allocate
   PM=repmat(NaN,1,12); % to hold gridpointmeans
   PS=repmat(NaN,1,12);  % to hold gridpoint std devs
   W5=repmat(NaN,1,ns1); % to hold weights on each station
   
   if inone==0; % if have at least one station in search radius
      
      % Rename unadjusted long-term means and sdevs
      TMN=UMN;
      TSD=USD;
      
      % Mean
      L1=~isnan(TMN'); % logical to non-NaNs of long-term means
      
      % Weight the long-term neans to gridpoint
      [WWM,NWM,JWM]=wgtdist1(A,L1,dcrit,dsrch,nmax,kk);
      % Note that all rows of WWM are identical. Only need first
      jcols=JWM(1):JWM(2);
      W5=WWM(1,jcols);
      wgtmean =W5;  % rename weights for the long-term means.  A rv.
      
      if strcmp(kdiag,'Yes');
         vlist=char(vlist,'wgtmean -- weights for interp of long-term means');
         eval(['save ' lookfile ' vlist wgtmean -append;']);
         clear wgtmean;
      end;
      
      L1=isnan(W5); % mark NaNs in W5
      sum5=sum(sum(L1)); % sum of NaNs in rv W5
      % Matrix mult with NaNs gives Nans, so change NaNs to zero
      if sum5>0;
         W5(L1)=0;
      end
      
      % Make temporary matrices TTMN and TTSD by replacing NaNs with zeros
      % in TMN, TSD before matrix mult.  Otherwise, NaNs in matrix mult, which
      % means all products NaN
      TTMN=TMN;
      TTSD=TSD;
      L1=isnan(TMN);
      sum5=sum(sum(L1));
      if sum5>0;
         TTMN(L1)=0;
      end
      L1=isnan(TSD);
      sum5=sum(sum(L1));
      if sum5>0;
         TTSD(L1)=0;
      end;
      
      PM=W5*TTMN;
      PS=W5*TTSD;
      clear WWM NWM JWM L1 sum5 TTMN TTSD
      
      %************* CONVERT GRIDPOINT STDZD ANOMALIES TO ORIGINAL UNITS
      
      disp('CALLING GRD2REG.M TO WEIGHT GRIDPOINT DATA TO REGIONAL DATA');
      IR2=[];
      [P,RA,R,RM,RS,JP3,wt2]=grd2reg(IR2,CG(nnn,:),Q,PM,PS,JP2);
      
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
   else; % no stns in search radius
      P=repmat(NaN,nrows,12);
   end; % if inone==0
   
      
   %******************** *****************************
   disp('BUILDING OUTPUT FILES')
   
   %****** OUTPUT THE RESULTS **************************
   % This section written and tested outside of regcli2.out, then
   % added later. 
   % 
   % RECALL THAT
   % B1 station filename list
   % C2 long lats of stations
   % CG long lats of gridpoints
   % NSTNS number of good stations matrix, all regions combined
   % P tsm of monthly gridpoint climate data
   % Q tsm ... as standardized anomalies
   % PM long term monthly means of gridpoint data
   % PS long term standard devs
   % TMN long-term (ref period) monthly means of station data
   % TSD long-term   std devs
   % W5 matrix of weights used to compute long-term gridpoint
   %	  means from long term station means
   % yr2 year vector for huge matrix
   % yr3 year vector for NSTNS after trimming off all-nan rows
   % filecode -- 1-char letter added to file prefixes for record keeping
   
   if inone==0; % if any stations in search radius
      % Reduce row size of P,Q by lopping off leading and trailing all-NaN years
      [mQ,nQ]=size(Q);
      i1=(1:mQ)'; % sequential index same row size as Q
      L1=isnan(Q);
      L1=(all(L1'))'; % cv, 1 if year has 12 NaNs, 0 otherwise
      i1(L1)=NaN;
      i1=trailnan(i1);
      i1=flipud(trailnan(flipud(i1)));
      if any(isnan(i1));
         i1(isnan(i1))=[];
      end;
      P=P(min(i1):max(i1),:);
      Q=Q(min(i1):max(i1),:);
      yr4=yr2(min(i1):max(i1));
      nyrs4=length(yr4);
   else; % no stations in search radius
      Q=[];
      P=[];
   end; % if inone==0
   
   
   % Gmn?.txt gridpoint interpolated monthly means
   % Gsd?.txt gridpoint interpolated std devs
   % Glist?.txt -- point #, lon, lat (5 station seq nos)-- these for means and std devs
   % Gwt?.txt point#, lon,lat, (5 distances (km)), (5 weights)-- these on means and std devs
   
   % ? == p, d, e, or m, depending on the type of climate data
   
   % Build gridpoint label
   Gname = [filecode int2str(nnn)];
   nombre=sprintf('%7s',Gname);
   nombre=strjust(nombre,'left');
   
   
   % Compute number of weights used for long-term mean interpolatoion
   if inone==0;
      w5=W5; % rv, weights 
      % find the non-zero elements
      j=find(w5~=0);
      w6=w5(j);
      nwgts=length(w6); % # weights
   else;
      w5=NaN;
      j=NaN;
      w6=NaN;
      nwgts=0;
   end;
   
   
   
   %  GLIST?.TXT  -- LIST OF GRIDPOINTS, LON,LAT, AND 'PREDICTOR' STATIONS
   if strcmp(kdiag,'Yes');
      fid1=fopen(['Glist0' ctype '.txt'],'w');
   else;
      fid1=fopen(['Glist' ctype '.txt'],'a');
   end;
   
   % Point list, with station filenames, long, lat, 
   fmt1='%7s %7.2f %6.2f  ';
      
   if strcmp(kdiag,'Yes');
      if nwgts>0;
         Btemp=B1(j,:); % names of climate station files used for weighting long-term means
      else;
         Btemp=NaN;
      end;
      for k=1:nwgts;
         if nwgts>0;
            btemp=Btemp(k,:);
            % Remove leading and trailing blanks
            btemp=strtok(btemp);
            btemp=fliplr(strtok(fliplr(btemp))); 
            Bnames{k}=btemp;
         else;
            Bnames{k}=NaN;
         end;
         
      end;
      vlist=char(vlist,'nwgts -- number of stations for intrp long-term means');
      vlist=char(vlist,'Bnames -- station filenames for long-term mean interp');
      eval(['save ' lookfile ' vlist nwgts Bnames -append;']);
      clear Bnames btemp Btemp;
   end;

   
   fprintf(fid1,fmt1,nombre,CG(nnn,1),CG(nnn,2));
   if nwgts>0;
      for k = 1:nwgts;
         str1=sprintf('%8s',B1(j(k),:));
         if k==nwgts;
            fprintf(fid1,'%s\n',str1);
         else;
            fprintf(fid1,'%s',str1);
         end;
      end
   else;
      fprintf(fid1,'%s\n',blanks(5));
   end;
   
   fclose(fid1);
   
   
   % GINDX -- CROSS-REFERENCE TO STN INDEX AND DISTANCE
   
   if strcmp(kdiag,'Yes');
      fid1=fopen(['Gindx0' ctype '.txt'],'w');
   else;
      fid1=fopen(['Gindx' ctype '.txt'],'a');
   end;
   fmt1='%7s ';
   fprintf(fid1,fmt1,nombre);
   if nwgts>0;
      for k = 1:nwgts;
         jindex = j(k);
         adist=A(jindex); 
         str1=sprintf('%4d ',jindex);
         str2=sprintf('(%5.1fkm)',adist);
         strall=[str1 str2];
         if k==nwgts;
            fprintf(fid1,'%s\n',strall);
         else;
            fprintf(fid1,'%s',strall);
         end
      end;
   else;
      fprintf(fid1,'%s\n',blanks(8));
   end;
   clear adist jindex str1 str2 ;
   fclose(fid1);

   
   
   
   %  GWGT.TXT --  WEIGHTS FOR PREDICTOR STATIONS 
    if strcmp(kdiag,'Yes');
      fid1=fopen(['Gwgt0' ctype '.txt'],'w');
   else;
      fid1=fopen(['Gwgt' ctype '.txt'],'a');
   end;
   fmt1='%7s ';
   fprintf(fid1,fmt1,nombre);
   if nwgts>0;
      for k = 1:nwgts;
         str1=sprintf('%7.4f ',w6(k));
         if k==nwgts;
            fprintf(fid1,'%s\n',str1);
         else;
            fprintf(fid1,'%s',str1);
         end
      end;
   else;
      fprintf(fid1,'%s\n',blanks(8));
   end;
   fclose(fid1);
   
   
   % Gmn -- interpolated gridpoint means
   if strcmp(kdiag,'Yes');
      fid1=fopen(['Gmn0' ctype '.txt'],'w');
   else;
      fid1=fopen(['Gmn' ctype '.txt'],'a');
   end;
   fmt1='%7s ';
   fmt2=' %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n';
   fprintf(fid1,fmt1,nombre);
   if nwgts>0;
      fprintf(fid1,fmt2,PM');
   else;
      fprintf(fid1,'%s\n',blanks(8));
   end;
   fclose(fid1);
   
   % Gsd -- interpolated gridpoint std dev
   if strcmp(kdiag,'Yes');
      fid1=fopen(['Gsd0' ctype '.txt'],'w');
   else;
      fid1=fopen(['Gsd' ctype '.txt'],'a');
   end;
   fmt1='%7s ';
   fmt2=' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n';
   fprintf(fid1,fmt1,nombre);
   if nwgts>0;
      fprintf(fid1,fmt2,PS');
   else;
      fprintf(fid1,'%s\n',blanks(8));
   end;
   fclose(fid1);
   
   
   % OUTPUT .mat files -- only if not diagnostics mode
   flname=[path5 Gname];
   if strcmp(kdiag,'Yes');
      vlist=char(vlist,'flname -- path\filename of output .mat point or grid series');
      eval(['save ' lookfile ' vlist flname -append;']);
   else;
      if inone==0;
         Z=[yr4 P];
         eval(['save ' flname  ' Z;']);
      else;
         Z=[];
         txt='No stations were in the search radius for this point';
         eval(['save ' flname  ' Z txt;']);
      end;
      
      % Set the filename prefix in Gnm in siteinf.mat, siteinfd.mat or whatever
      eval(['load ' sitefile ' Gnm;']);
      Gnm{nnn}=Gname;
      eval(['save ' sitefile ' Gnm -append;']);
   end; % if strcmp(kdia,'Yes')
   
end; % for nnn= loop over gridpoints


seconds=etime(clock,t0);
tminutes=seconds/60;
clc;
disp('Finished running regcli5.m');
disp(['  Total of ' num2str(npnts) ' stations']);
disp(['  Elapsed time = ' num2str(tminutes)]);

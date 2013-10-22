function regcli3a
% regcli3a:  interpolation of station data to points by standardized anomaly
% CALL: regcli3a
%
% Meko 11-21-97
%  Revised 7-24-98: zero as interpol pcp if all predictor stations pcp zero
%
% Application.  Need to get monthly files of tmp at tree locations or
% gridpoints from GHCN station data.  Calls regcli3b.m to do most of work
% 
%
%************  IN ******************************
% User prompted for whether this is a "massive urban heat island T" or "other" run
%   If "other", user prompted for whether data pcp or tmp.  Response will determine
%   naming of output .mat data files.  For the "massive" run, you get files numbered
%   TD1.mat, TD2.mat, ...
%   For the "other" runs, you get .mat files numbered
%    pnt1t.mat, pnt2t.mat,...,  or pnt1p.mat, pntp.mat, ...
%   
% User prompted for name of .mat file with following contents
%
% datin{}
%
%	{1} paths {1 x 4} paths to relevenat files
%		1 - working directory (e.g., c:\projs\wrk0\)
%		2 - path to put output files in (e.g., c:\projs\wrk0\outfls\)
%		3 - path to the dir with climate database (or wherever you have
%         stored lat-lon files, elev files, ...
%			(e.g., c:\data\foxy\ghcnt\tmaxadj\)
%		4 - path to dir with individual station input .mat monthly clim files
%			(e.g., c:\data\ghcn\tmpv2_0\matflsd\
%	{2} files {1 x 4} filenames for:
%			1 - ascii information file for points (in path 1). This info
%				originally formatted to accomodate production of the file from
%				my foxpro tree-ring "beef.dbf" database.  For other applications,
%				such as interpolating to a lat-long point, can supply surrogate
%				info in the required character slots.  Column:
%
%           3-10  name of .mat file holding tree-ring chron; if
%				16-35 truncated point or site name
%				47-48 state code
%				53-55 country coude
%				61-64 species code
%				69-72 el (m)
%				
%			2 - lon-lat file for points (in path 1).  Two values per line. Long
%				in col 1, lat in col 2.  Format is decimal degrees, with W long neg
%
%			3 - ascii info file for stations (in path3)
%				Format of this file originally written to accomodate automatic
%				building from foxpro ghcn version 2 temperature database.If
%				data from other source, manually build to match following
% 				format:
%				
%				col 3-5 country code (e.g., USA, MEX, CAN)
%				col 11-30 first 20 chars of station name. MAKE SURE you make
%               lines at least 30 chars wide
%
%			4 - ascii file of elevations (2-col) for stations (in path 3).
%				Two cols to accomodate automatic building from the version 2 
%				ghcn temperature database.  That database has a recorded and
%				interpolated elev.  Program will use the recorded value unless
%				-999, when it will use the interpolated.  For other data 
%				sources, just put the elev in col 1 and a dummy col or -999
%				in col 2.  El in meters.
%
%			5 - lon-lat file for stations (in path 3). Format
%				same as for lon-lat file for points
%
%			6 - ascii file with names of station .mat files (in path3) 
%				One filename (e.g. tucson) per line.  No spaces inside the name.
%				No need for trailing .mat.   Name need not begin in col 1.
%				
%	{3} (2 x 2)i start, end year for 
%		row 1: reference period
%		row 2: big matrix period
%	{4} (1 x 3)r  distance info
%			dcrit -- critical threshold distance to avoid over-weighting site in
%				distance weighting (km)
%			dsrch -- search distance for weighting (km)
%			nmax --  max number of stations to weight to a gridpoint
%	{5} (1 x 4)i   options
%			kopts(1): mode for running program
%					==1: diagnostic (lots of feedback)
%					==2: operational (just the facts, Ma'am)
%			kopts(2):==1 --> subst. zero for negative monthly regional values
%				    	==2 --> accept negative monthly regional values
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
%
%	{6} thresholds
%			pgood (1 x 1)r  minimuma acceptable decimal proportion of monthly data
%					not NaN for station to be considered as candidate for use in
%					interpolation (typically, .75).  This proportion is in the
%					reference period.
%
%
%****************** OUT ********************************************
%
% No output args.  Just calls regcli3b.m, which generates output for each point
%
%
%******************* STEPS ********************************************
%
% Read input and unpack cell info
% Get number of points, stations from  from row sizes of lon-lat files
% Compute distance matrix, points to stations
% Get a single elevation value for the stations
% Read ascii text info for stations; country and station name; condense info and
%      store in string matrix
% Read some ascii text info for points, ...
% Loop over stations, storing .mat file prefixes and
%		building logical pointer to stations with adequate time coverage in reference period
% Concatenate station info matrix, including long, lat, elev, year range, etc
% Build matrix of path/filenames to input .mat climate files
% Cull submatrices based on screening for adequate time coverage in ref period:
%		-rows of station lon-lat 
%		-rows of station info file
%		-rows of string matrix with .mat file names for stations
%		-cols of point-to-station distance matrix
% Loop over points, calling regcli3b.m to get  output files
%    (See comments in regcli3b.m for technical details of the stdzd anomaly method)
%
%************* USER FUNCTIONS CALLED
% 
%
% REGCLI3B M          26,625  06-21-98  3:40p Regcli3b.m
% GCDIST   M           2,623  07-23-96 11:07a GCDIST.M
% DEC2DMS2 M           2,073  04-30-96  5:04p DEC2DMS2.M
% FPFMON1  M* not needed, but nice to have to format output
% DMS2RAD  M             284  09-08-93 11:15p DMS2RAD.M
% JUXTA1   M           1,285  09-09-93  1:46p JUXTA1.M
% KEYDIST  M             734  06-05-98  2:10p KEYDIST.M
% STDNAN   M           2,700  01-21-97  9:03a STDNAN.M
% MEANNAN  M           1,477  01-21-97  8:42a MEANNAN.M
% WGTDIST2 M           5,556  11-24-97  1:25p WGTDIST2.M
% NONAN1   M           1,803  11-14-97 10:07p NONAN1.M
% SUMNAN   M             940  01-28-97  4:55p SUMNAN.M
% TRAILNAN M             649  02-28-97  4:38p TRAILNAN.M

%********* Determine data type
kmassive = questdlg('Is this a massive run to get all version 2 urban t files?');
if strcmp(kmassive,'Yes');
   preftype = []; 
else
   kother= menu('Data Type','Precipitation','Temperature');
   if kother==1;
      preftype='p';
   else
      preftype='t';
   end
end

%********* UNPACK CELL INPUT AND CHECK INPUTS

% Get name of .mat file with input cell
[filetemp,pathtemp]=uigetfile('*.mat','Input .mat file with datin{}, etc');
pftemp=[pathtemp filetemp];
eval(['load ' pftemp]); 
clear filetemp pathtemp pftemp
% Now have the cell info datin{} in workspace

paths=datin{1};
files=datin{2};
pdref = datin{3}(1,1:2);
pdstore=datin{3}(2,1:2);
dcrit=datin{4}(1);
dsrch = datin{4}(2);
nmax = datin{4}(3);
kopts=datin{5};
pgood=datin{6};


%*************** STORE LON/LAT MATRICES, AND COMPUTE NUMBER OF STATIONS AND POINTS

%------------ STATIONS
pftemp=[paths{3} files{5}];
eval(['load ' pftemp ';']);
eval(['xystn = ' strtok(files{5},'.') ';']);
[m1,n1]=size(xystn);
clear pftemp
% xystn has the lon-lats for stations
% m1 is the number of stations


%------------ POINTS
pftemp=[paths{1} files{2}];
eval(['load ' pftemp ';']);
eval(['xypnt = ' strtok(files{2},'.') ';']);
[m2,n2]=size(xypnt);
clear pftemp
% xypnt has the lon-lats for points
% m2 is the number of points

%--------------  COMPUTE POINT-TO-STATION DISTANCE MATRIX (KM)
A = gcdist(xypnt,xystn);
[mA,nA]=size(A);

% ------------- GET A SINGLE ELEV VALUE (M) FOR EACH STATION
pftemp=[paths{3} files{4}];
eval(['load ' pftemp ';']);
eval(['elev1 = ' strtok(files{4},'.') ';']);
[mtemp,ntemp]=size(elev1);
clear pftemp
% elevstn has the lon-lats for stations; first col is given el, second
% is and interpolated value.  Want to accept the given, but if the 
% given is -999, to use the interpolated
elevstn=elev1(:,1);
Ltemp=elev1(:,1)==-999;
if any(Ltemp);
   elevstn(Ltemp) = elev1(Ltemp,2);
end
clear elev1 Lemp mtemp ntemp pftemp
% UPDATE: station elevations now in elevstn



%**************  GET AND STORE SOME TEXT INFO FOR STATIONS
SS1= repmat(blanks(24),m1,1); % to store country and station name


pftemp=[paths{3} files{3}];
fidtemp = fopen(pftemp,'r');
for m = 1:m1; % loop over stations
   c=fgetl(fidtemp);
   c1 = c(3:5); % country code
   if strcmp(c1,'403') | strcmp(c1,'CAN');
      ctry = 'CAN';
   elseif strcmp(c1,'425') | strcmp(c1,'USA');
      ctry = 'USA';
   else;
      fclose all
      error('So far, only USA and Canada country code recognized');
   end
   c2 = c(11:30); % first 20 chars of station name
   SS1(m,:)=[c2 ' ' ctry];
end
fclose(fidtemp);
clear fidtemp c c1 c2 ctry
% UPDATE: SS1 has country and station name 


%**************  GET AND STORE SOME TEXT INFO FOR POINTS
SP1= repmat(blanks(80),m2,1); % to store country and station name


pftemp=[paths{1} files{1}];
fidtemp = fopen(pftemp,'r');
for m = 1:m2; % loop over points
   c=fgetl(fidtemp);
   c0 = sprintf('%4.0f',m);
   c1 = c(3:10); % name of .mat file
   cname = blanks(8);
   c1 = strtok(c1,'.');
   len1 = length(c1);
   cname(1:len1)=c1;
   c2 = c(16:35); % truncated site name
   c3 = c(47:48); % state
   c4 = c(53:55); % country
   c5 = c(61:64); % species
   c6 = [c(69:72) 'm']; % elev (m)
   clon = sprintf('%7.2f',xypnt(m,1));
   clat = sprintf('%5.2f',xypnt(m,2));
   str1 = [c0 '-' cname ' ' c2 ' ' c3 ' ' c4 ' ' c5 ' '];
   str2 = [clon ' ' clat ' ' c6];
   str3 = [str1 str2];
   nlen = length(str3);
   SP1(m,1:nlen)=str3;
end
fclose all;
clear c c0 c1 c2 c c4 c5 c6 clon clat str1 str2 str3 nlen cname m
% UPDATE:  SP1 has text info for the m2 points



%********** LOOP OVER STATIONS, STORING .MAT FILE NAMES, AND BUILDING LOGICAL
% POINTER TO STATIONS WITH ADEQUATE COVERAGE FOR REFERENCE PERIOD

% Read and store .mat station filenames
SS2 = repmat(blanks(8),m1,1);
YRS = zeros(m1,2); % to hold first, last year with any monthly data

pftemp=[paths{3} files{6}];
fidtemp = fopen(pftemp,'r');
for n = 1:m1;
   c = fgetl(fidtemp);
   c=deblank(c);
   c=deblank(fliplr(c));
   c=fliplr(c);
   clen = length(c);
   SS2(n,1:clen)=c;
end

disp('Checking time coverage in .mat files');
yrref = (pdref(1):pdref(2))';
nref = length(yrref);

Luse = logical(zeros(m2,m1)); % m2 points by m1 stations

for n = 1:m1; % loop over stations (cols of Luse)
   flnm = [paths{4} strtok(SS2(n,:))];
   eval(['load ' flnm]);
   Zdata = X(:,2:13);
   yrZ = X(:,1);
   YRS(n,:)=[min(yrZ) max(yrZ)];
   L1 = yrZ>=pdref(1) & yrZ<=pdref(2);
   if sum(L1)>0.75*nref;
      Luse(:,n)=1;
   end
end
Luse=logical(Luse); % Note that rows of Luse are all the same
clear Zdata yrZ L1 flnm clen


%**************  CONCATENATE  STATION INFO MATRIX
SS3=repmat(blanks(80),m1,1);


for n=1:m1;
   str1 = sprintf('%7.2f %5.2f ',xystn(n,:)); % lon and lat
   str2 = sprintf('%4.0fm ',elevstn(n)); % elev (m)
   str3 = sprintf('%4.0f %4.0f ',YRS(n,:));
   strold = [SS1(n,:) ' ' SS2(n,:) ' '];
   str4 = [strold str1 ' ' str2 str3];
   len1 = length(str4);
   SS3(n,1:len1)=str4;
end
clear str1 str2 str3 strold str4 len1

   
%*********************  BUILD MATRIX OF PATH/FILENAME TO CLIMATE FILES
Sfiles = repmat(blanks(70),m1,1);
nbig = 0;
for n = 1:m1;
   pftemp = [paths{4} SS2(n,:)];
   nlength = length(pftemp);
   nbig = max(nlength,nbig);
   Sfiles(n,1:nlength)=pftemp;
end
Sfiles=Sfiles(:,1:nbig);
clear nbig nlength pftemp n   



%********************* CULL SUB-ROWS OF CLIMATE  AND RELATED MATRICES 
% Recall that already have Luse pointing to stations that have required
%   pctg of non-missing data in the reference period
% Get submatrices with rows or cols for stations with insufficient
% time coverage in the refrence period removed

% Now Figure that will not want to pass any more than 20 stations as 
% potential predictors.  And if kopts(3)==2, will only need to consider
% stations within dsrch km of the point.  Build logical pointer matrices to
% the acceptable stations
dkey = keydist(A,min(m1,20)); % cv of distance to 20th nearest station for each point
%  or to the m1 th nearest if total number of stations less than 20
DKEY = repmat(dkey,1,nA); % dupe to matrix with col size equal to number of stations
Lkey = A<=DKEY; % 1 if distance within dist to 20th nearest station
Lsrch = A<=dsrch; % 1 if distance less than search radius

% Build logical col pointer
if kopts(3)==1; % not restricted to stations in dsrch
   LOK =  Luse & Lkey;
else; % kopts(3)==2, so stations must be in search radius
   LOK =  Luse & Lkey & Lsrch;
end




%******** CALL REGCLI3B.M POINT BY POINT TO COMPUTE THE INTERPOLATED VALUES
disp('Making Calls to To regcli3b.m for each point');

for n = 1:m2; % loop over points
   disp(['   Point # ' int2str(n)]);
   miscin{1}=xypnt(n,:); % lon/lat of point
   
   if n ==34;
      disp('At site 34');
   end
   
      
   % Cull matrices 
   flsin=Sfiles(LOK(n,:),:); % input filenames
   miscin{2}=xystn(LOK(n,:),:); % lon/lat of stations
   miscin{3}= A(n,LOK(n,:)); % distance vector
   Sinf={SP1(n,:),SS3(LOK(n,:),:)};

   
   miscin{4}=[pdref;pdstore];
   miscin{5} = datin{5};
   miscin{6} = datin{4};
   miscin{7}={paths{4},paths{2}};
   miscin{8}=kmassive;
   miscin{9}=preftype;
   pntno = n;
   regcli3b(flsin,miscin,Sinf,pntno);
end

disp('OK')

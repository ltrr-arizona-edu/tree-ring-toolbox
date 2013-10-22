function  ptest(mss,P,kmode,t1,t2,Nd,Ndmask,A4)  
%
% Estimation of missing monthly precipitation or temperature data
%
% D. Meko,  10-6-94;  last revised, 7-28-96
%
% Model for ptest.m was zest.m.  Main differences:  (1) ptest.m
% has different input requirements, (2) ptest.m produces different 
% output, (3) ptest has more quality control, (4) ptest.m accepts
% NaN as missing values
%
%*********  INPUT -- all stored in a .mat prepared beforehand  ********
%
% Nd (1 x 1)i    number of input series
% mss  (1 x 1)r  missing value code [say, -9.99]
%		NaN is also recognized as missing value; if NaN is used, 
%		whatever	you specify for mss will be checked for in addition
%		to NaN's
% P    (1 x 1)L   either precip (1) or temperature (0)
% t1  (1 x 1)i   earliest beginning year of any series
% t2  (1 x 1)i  latest ending year of any series 	
% kmode (1 x 1)i  1== srp mode, Nd, Ndmask, A4 not passed as args, 
%		but built from the ascii info file 
% Ndmask (1 x Nd)L  1=estimate for this series;  0=skip
% A4 (Nd x 10)i
% 	sequence numbers of predictors, space separated, ranked in 
%		priority left to right; zero-filled
%
% mat files with input 13-col monthly ppt data.  One file for each
% of Nd stations.  No headers.  Missing values as mss or NaN. These mat
% files must be available in the path or current directory. For Meko 
% system, root name is 3-5 characters, followed by "P" or "T".  Thus
% ALPINP.MAT for Alpine,  BISB2P.MAT for Bisbee 2, AJOP.MAT for Ajo.
% A matrix "Z" in the .mat file is assumed to hold the pcp data.
%
%
%************** ADDITIONAL INPUT FROM AN ASCII FILE ****************
%
% 
%***************   OUTPUT ****************************
%
% Two new .mat files of "data" are produced
% for each input series.  File with "X" appended to name is the 
% original data with estimates substituted for missing values.
% File with "E" appended to name holds sequence numbers of 
% predictors for each predicted value.  Thus, might have
% AJOP.MAT -->  AJOPX.MAT, AJOPE.MAT.
% Special meaning is given to 0 0r -1 in the "E" file.  0 means no
% no estimation done.  -1 means "zero-type" estimate (see notes in 
% meanrt2.m
%
% Files with the root name of the station, ending in "M", "S", and
% "R" are produced.  These files hold the mean ratios or mean
% departures, the sample sizes for computing the mean ratios or
% departures, and the correlation coefficients between the
% monthly series of the key station and each predictor.
%
% If no data need to be estimated for the key station, or if
% the station has been masked (its Ndmask entry is 0), then no
% output data files are produced.
%
%
%*********************** NOTES **********************************
%
%*** PROGRAM STEPS
%- Initialize a matrix X1 covering years (rows) t1 to t2 and having
%  12 columns for each station.  A column holds the monthly pcp for
%  a station
%
%
%***************  USER-WRITTEN FUNCTIONS NEEDED 
%
% meanrt1.m    calc mean ratios for predictors/key
% meanrt2.m    apply the mean ratios to estimate vector of
%	missing values for a specific month of year at one
%	key station.
%
%*******************  END OF COMMENTS  ***************************

%****************************************************************
%

% Normally, if no data need to be estimated for the key station,
% no need to make a files with mean ratios, mean departures, and
% correlations between monthly series of key and predictor stations.
% But for some exploratory studies, you may simply want to find
% out how well a station's monthly pcp or tmp records relate to 
% those at some nearby stations.  So an option exists to make
% the "M", "S" and "R" files even when no data need to be
% estimated at the key station
disp('You might want M, S and R files even if no data needed ')
disp('to be estimated at the key station.  If so, answer "Y" to')
disp('the following prompt')
s1 = input ('M,S and R files anyway (Y/[N])?','s');
if isempty(s1) | s1=='n',
	s1=='N';
elseif  (s1=='Y' | s1=='y')
	s1='Y';
else
	error('Invalid response to s1 prompt')
end

clc

%**********************************************************
% Get the file names of monthly pcp .mat files. Store the names
% in matrix A2, 8 cols, with name left justified
blanky='        '; % blank vector
A2=blanky(ones(100,1),:); % Blank matrix, initialized
[file1,path1]=uigetfile('*.txt','filename file -- e.g., nearaz.txt');
pf1=[path1 file1];
fid1=fopen(pf1,'r');
% file structure
%    cols 1-2, right-justified seq number;  key station is # 1
%    cols 4-11, filename, without suffix, might be left or right just
%    cols 13-32, station name, left just
%    cols 34-40, distance (km) key to this, right just
%	  cols 42-45 first year of monthly data
%	  cols 47-50 last year of monthly

if kmode==1; % srp mode
	k7=1;
	ncount=0;
	while k7==1;
		c=fgetl(fid1);
		ncount=ncount+1;
		n=ncount;
		if feof(fid1);
			k7=0;
			Nd=ncount-1;
		else
			aname=c(4:11);
			% get rid of leading zeros, if any
			aname=fliplr(aname);
			aname=deblank(aname);
			aname=fliplr(aname);
			% store in A2
			nchar=length(aname);
			A2(n,1:nchar)=aname;
			% want earliest start, latest end year of series
			tgo = str2num(c(42:45));
			tsp = str2num(c(47:50));
			if ncount==1;
				t1=tgo;
				t2=tsp;
			else
				t1=min([t1 tgo]);
				t2=max([t2 tsp]);
			end
		end
	end
	A2=A2(1:Nd,:);
	Ndmask=zeros(1,Nd);
	Ndmask(1)=1;
	A4=zeros(Nd);
	A4(1,:)=[2:Nd 0];
end


fclose(fid1);


% STORE MONTHLY DATA FOR ALL STATIONS IN ONE BIG TSM;  FIRST STATION'S
% DATA IN COLS 1-12, SECOND'S IN COLS 13-24, ETC

%Initialize the storage matrix X1
a=NaN;
mX1=t2-t1+1;  % Number of rows
nX1=Nd*12;   % Number of cols
X1=a(ones(mX1,1),ones(nX1,1));  % initialize X1 as NaN

% Initialize storage for first, last years of monthly series
T=a(ones(Nd,1),ones(2,1));




% Get monthly data and put in big matrix
for n = 1:Nd;  % loop over stations
	% get the root of filename
	fn=A2(n,:);  % Root part of name, without suffix .mat
	fn=strtok(fn,' ');

	% load file of monthly data; put data in W, year in yr
	eval(['load ' fn]);
	eval(['W = Z;']);  % store the series in W
	yr = W(:,1);   % year cv
	txt1=['Series ',int2str(n),'  ',fn]; % will use in error messages

	% start and end year of monthly data for this station
	t3=min(yr);
	t4=max(yr);

	% Monthly data should have been prepared so that no entire
	% years (rows) are missing from the monthly data.  Check 
	% anyway to see that first and last years in col 1 are 
	% consistent with the number of rows in the monthly data
	% Check that start, end years are inclosed by years t1,t2
	% specified for master matrix.  Check also for 
	% other inconsistencies in years.

	% No station's monthly series should be outside period
	% specified for master storage matrix
	if t3<t1 | t4>t2,
		disp(txt1)
		error('Years outside range specified by t1,t2')
	end
	

	if size(W,1) ~= (max(yr)-min(yr)+1);
		disp(txt1);
		error('Row size inconsistent with start, end years')
	end

	% Year vector should be monotonic, increasing by one
	dyr = diff(yr);
	if ~all(dyr);
		disp(txt);
		error('Year vector not monotonic increasing by one')
	end

	% OK, years info passed the checks.
	% Store start, end years ;  will need these later
	T(n,1)=t3;
	T(n,2)=t4;

	% start, end cols of X1 for this stations monthly data
	j1=n*12-11;
	j2=j1+11;

	% rows in X1 for this station's data
	iX1= yr-t1+1;   % 

	% Store the monthly data in X1
	X1(iX1,j1:j2) = W(:,2:13);

end;   % of n= loop over stations
disp(['Monthly data for the ' int2str(Nd) ' stations has ']);
disp('been stored in the matrix X1');
disp('');
disp(['Year Range of X1 is ' int2str(t1) '-' int2str(t2)]);
disp('Next will start loop over key stations');
disp(' ')
disp(' ')


% Loop over key station to estimate data

for n=1:Nd
	if Ndmask(n);  % no need to estimate for masked stations
	fn=A2(n,:);  % Root part of name, without suffix .mat
	fn=strtok(fn,' ');
	disp(['Starting estimation for station ',int2str(n),' ',fn]);
	txt1=['Series ',int2str(n),'  ',fn]; % will use in error messages

	% Start and end years of monthly data for this series
	t3=T(n,1);
	t4=T(n,2);

	nyrs = t4-t3+1; % number of years of monthly data for this station
	yr = (t3:t4)'; % year vector for station

	%  Allocate storage for: 
	D = ones(nyrs,12);  % data block, same size as monthly data
	%	monthly data
	Y1 = [yr D];   % for monthly data with estimates filled in
	Y2 = [yr zeros(size(D))]; % for estimator-station pointer
	M=zeros(10,13); % mean ratios or departures for up to 10 predictor
	S=zeros(10,13); % sample sizes for mean-ratio or mean departure
	R=zeros(10,13); % correl coefs monthly series for key and preds

	% Make pointer to cols of X1 with data for key stn and predictors
	% Pull out block of data covering the key-station years
	jk = n;  % sequence number for key station
	jp = A4(n,:);  % sequence nos for predictor stations;
	%  right-filled with zeros
	jp = jp(jp~=0);  % predictor sequence numbers
	np = length (jp);  % number of predictor stations
	jall = [jk jp]; % seq numbers of key station and predictors


	% Loop over the 12 months of the year for this key station
	for j = 1:12;
		jcols = (jall*12)+j-12 ;  % cols of desired data in X1
		iX1= t3-t1+1;   % row-index for start of segment to pull
		iX2= t4-t1+1;   % row-index for end of segmnt
		B1 = X1(iX1:iX2,jcols); % data for key and predictor stations
		% for period of key station

		% Compute mean ratios; mean departures; the number of valid
		% years the ratios or departures are based on, and the 
		% correlations between monthly data for key and predictor
		% stations; recall that np is the number of predictors,
		% mss is the missing-value code, and P is a logical scalar --
		% 1 for pcp, 0 for tmp.
		[mr,md,nomr,ss,rr]=meanrt1(B1,np,mss,P);
		

		% Matrix M will hold the mean ratios or departures,
		% S will hold the sample sizes the ratios or departures are
		%   based on,
		% R will hold the correlation coeffficients measuring
		%   strength of relationship of monthly series

		% Put predictor sequence number in column 1 of M,S,R
		M(1:np,1) = jp';
		S(1:np,1) = jp';
		R(1:np,1) = jp';

		% meanrt1.m returned both the mean ratios and mean departures.
		% Use the appropriate one depending on the data type
		% (pcp or tmp)
		if P; % If datat-type is pcp use mean ratio
			M(1:np,j+1) = mr';
		else; % if data-type use mean departures
			M(1:np,j+1) = md';
		end

		% Store sample sizes and correlation coefficients
		S(1:np,j+1) = ss';  % sample size
		R(1:np,j+1) = rr'; % correlation coef

		% Estimate missing monthly values using mr or md
		[y1,y2] = meanrt2(B1,jp,mr,md,nomr,mss,P);

		Y1(:,j+1) = y1; % Monthly data, original if exists, with missing
			% values replaced by estimates 
		Y2(:,j+1) = y2; %  Vector same size as data, telling which
			% predictor used for that value;  zero if no estimate needed
	end;  % of for n= loop over stations

	% Trim the excess of the 10 rows from the files holding the
	% ratios or departures, sample sizes, and correlation coeffs
	M=M(1:np,:);
	S=S(1:np,:);
	R=R(1:np,:);

	% save the files with the station's estimated data,
	% the predictor identificaton, mean ratios/departures,
	% and the sample sizes for mean ratios/departures

	% Append appropriate letter to file name for output type
	fnoutx = [fn 'X';]; % The monthly data with estimates substitued
		% for missing values
	fnoute = [fn 'E']; % Sequence number for estimator stns
	fnoutm = [fn 'M']; % Mean ratios or mean departures, depending on P
	fnouts = [fn 'S']; % The sample size for mean ratios or departures
	fnoutr = [fn 'R']; % The correlation coefficients

	% Want an 'E' file for a staion only if any data were estimated.
	% The E file will contain only a subset of rows (years) -- 
	% corresponding to the years for which estimates are needed.
	% If no data needed to be estimated, also no need to make an
	% "X" file, "M", "S", or "R" file

	Q=Y2(:,2:13);  % predictor seq numbers, or 0 if no prediction,
		% or -1 if special-case zero-rainfall prediction
	L = (any(Q'))';
	if sum(L)~=0;  % at least one year has some estimated data
		Y2=Y2(L,:);
   	eval (['save ' fnoute  ' Y2']);
		eval (['save ' fnoutx   ' Y1']);
	else
		disp('No data needed to be estimated for above station')
	end

	if (sum(L)~=0) | (s1=='Y');
		eval (['save ' fnoutm  ' M']);
		eval (['save ' fnouts  ' S']);
		eval (['save ' fnoutr  ' R']);
	end

	else;  % if Ndmask to skip this series
	end; % if Ndmask to skip this series

end; %  'n=' loop over stations
	
disp('YA DONE GOOD, BUBBA!')

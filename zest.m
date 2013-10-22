% zest.m  
%
% Estimation of missing monthly precipitation or temperature data
%
% D. Meko,  10-6-94;  last revised, 10-25-94
%
%
%*********  INPUT -- all stored in a .mat prepared beforehand  ****************
%
% Nd (1 x 1)i    number of input series
% mss  (1 x 1)r  missing value code [say, -9.99]
% stage (1 x 1)i   estimation step;  will replace col 5 of input file
%	names.  Say an input deck is lem3x.mat, and stage=1; output
%	file s lem31.mat
% P (1 x 1)L   either precip (1) or temperature (0)
% t1 (1 x 1)i   earliest beginning year of any series
% t2 (1 x 1)i  latest ending year of any series 	
% Ndmask (1 x 50)L  1=estimate for this series;  0=skip
% A2 (Nd x 25)s      one row for each station;   vbls not named in file
%   1-3   sequence number of this station 
%   11-15  mat-file name (prefix) of input data series
% A3 (Nd x 23)i   one row for each station;  vbles not named in file
%   seq  numeric sequence number (matches that in cols 1-3 of A2;
%   t3  beginning year of this series
%   t4  ending year
%   ... seq numbers of predictor station, up to 20, space-separated
%
% mat files with input 13-col monthly ppt data.  One file for each
% of Nd stations.  No headers.  Missing values as mss.  These mat
% files must be available in the path or current directory
%
%
%
%***************   OUTPUT ****************************
%
% mat files corresponding to each input series.  Two files produced
% for each input series.  File with col 5 of name changed is the file
% with original data with estimates substituted for missing values.
% File with col 5 changed and with "E" appended to name is file with
% sequence numbers of predictors used for each value.  0 means no
% no estimation done.  -1 means "zero-type" estimate (see notes in 
% meanrt2.m
%
%
%
%***************  USER-WRITTEN FUNCTIONS NEEDED 
%
%  meanrt1.m    calc mean ratios for predictors/key
%  meanrt2.m    apply the mean ratios to estimate vector of
%	missing values for a specific month of year at one
%	key station.
%
%*******************  END OF COMMENTS  ***************************

mX1=t2-t1+1;  % Number of rows in combined data matrix
nX1=Nd*12;   % Number of cols
X1=mss(ones(mX1,1),ones(nX1,1));  % initialize X1 with missing values


% load files of station ppt or T data
for n = 1:Nd;  % loop over key stations

   fn=A2(n,11:15);  % filename of input series (a mat-file)
   t3=A3(n,2);  % start year of key station
   t4=A3(n,3);  % end year of key station

   %  start, end cols of X1 that series will go into
   j1=n*12-11;
   j2=j1+11;

   % load file containing the series
   eval(['load ' fn]);
	eval(['W = X;']);
   %eval(['W= '  fn ';']);   % store the series in W
   yr = W(:,1);   % year cv
   iX1= yr-t1+1;   % row-index for putting series in X1
   X1(iX1,j1:j2) = W(:,2:13);  % put the station's data in X1

end ;   % of for-loop putting station series into combined matrix X1
disp('All series have been placed into X1');
disp('Next will start loop over key stations');
pause(2);


% Loop over key station to estimate data

for n=1:Nd
   disp(['Starting estimation for station ',int2str(n)]);
   if Ndmask(n);
   fn=A2(n,11:15);  % input data file name
   t3=A3(n,2);  %start year
   t4=A3(n,3);  %end year
   nyrs = t4-t3+1;
   yr = (t3:t4)';

%  Allocate
   Z = ones(nyrs,12);  % allocate
   Y1 = [yr Z];   % will hold the estimates
   Y2 = [yr zeros(size(Z))]; % will hold seq number of estimator

   % Make pointer to cols of X1 holding data for key and predictor stations
   % Pull out block of data covering the key-station years
   jk = A3(n,1);  % sequence number for key station
   jp = A3(n,4:23);  % sequenc nos for predictor stations (zero-filled)
   jp = jp(jp~=0);  % get rid of fill-zeros
   jall = [jk jp];
   np = length (jp);  % number of predictor stations

   % Allocate
   M=zeros(10,13); % Will hold mean ratios or departures
   S=zeros(10,13); % will hold sample size that mean-ratio or mean-dep
   R=zeros(10,13); % will hold correl coefs
  		% based on
   % Loop over the 12 months of the year for this key station
   for j = 1:12;
      jcols = (jall*12)+j-12 ;  % cols of X1
      iX1= t3-t1+1;   % row-index for start of segment to pull
      iX2= t4-t1+1;   % row-index for end of segmnt
      % Pull the segment
      B1 = X1 (iX1:iX2,jcols);  

      % Compute mean ratios and mean departures
		 [mr,md,nomr,ss,rr]=meanrt1(B1,np,mss,P);
		
		% Store mean ratios and their sample size
		M(1:np,1) = jp';
		S(1:np,1) = jp';
		R(1:np,1) = jp';

		if P
			M(1:np,j+1) = mr';
		else
			M(1:np,j+1) = md';
		end
		S(1:np,j+1) = ss';  % sample size
		R(1:np,j+1) = rr'; % correlation coef

		% Estimate missing monthly values using mr or md
		[y1,y2] = meanrt2(B1,jp,mr,md,nomr,mss,P);

		Y1(:,j+1) = y1;
		Y2(:,j+1) = y2;
	end

	% save the files with the station's estimated data,
	% the predictor identificaton, mean ratios/departures,
        % and the sample sizes for mean ratios/departures
	fnout = [fn(1:4) int2str(stage)];
   	fnoute = [fnout 'E'];
	fnoutm = [fnout 'M'];
	fnouts = [fnout 'S'];
	fnoutr = [fnout 'R'];	

	eval (['save ' fnout  ' Y1']);
   	eval (['save ' fnoute  ' Y2']);
	eval (['save ' fnoutm  ' M']);
	eval (['save ' fnouts  ' S']);
	eval (['save ' fnoutr  ' R']);
   else
   end
end
	
disp('YA DONE GOOD, BUBBA!')

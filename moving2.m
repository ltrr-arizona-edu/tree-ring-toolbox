function D=moving2(filesin,jcol,YRS,kopt,wsize1,wsize2,textin,seasons,thresh,contr);
% moving1:  color map of moving sum, average, median, or frequency of events
% CALL: D=moving2(filesin,jcol,YRS,kopt,wsize1,wsize2,textin,seasons,thresh,contr);
%
% Meko 4-30-97
%
% Differs from moving1.m in that moving2.m attempts to use 'actual' data
% for startup of any windowed analysis period. If I get moving2.m working,
% want eventually for it to supercede moving1.m
%
%******************* IN *****************************************************
%
% filesin {}  cell array, 1 row, containing paths/names of files holding input
%   time series data.  col size of filesin tells how many files hold input data
% jcol (1 x ?)i  sequence number or index to cols of the files, depending on 
%   value of kopt(1), telling where to the data for the time series:  
%   If kopt(1)==1, meaning annual or single-season mode, a single time series is needed,
%      filesin should contain only one file name, and jcol should be a scalar telling
%      which col (in matrix X) the series is to come from.  
%   If kopt(1)==2, meaning multi-season mode,
%      time series must be built from individual season series, which each are 
%      assumed to come from a different file.  Entries of jcol say which column in
%      each file to get the data from.  For example, if have two files, and jcol is
%      [3 3], means get season 1 from col 3 in first file of filesin, and season 2 from
%      col 3 of second file of filesin.  
%   If kopt(1)==3, meaning monthly mode,
%      jcol is assumed to be a scalar indicating the sequence number of the 
%      desired series from possibly multiple groups of 12 monthly series in the 
%      input matrix. Thus jser==1 means use cols 2-13, while jcol=2 says use
%      cols 14-25, etc.
% YRS (2 x 2)i  Start and end years of 
%   row 1: period for plot and tabular summary of moving window quantities
%   row 2: period to compute means  (anomalies might be needed)
%   row 3: period for analysis (moving averages, etc., computed for this
%			period before zeroing in on summary and plot period indicated by 
%			YRS (1,:).
% kopt (1 x 6)i options
%   kopt(1) series type for analysis (see notes)
%        ==1 annual or single-season mode
%        ==2 multi-season mode
%        ==3 monthly mode
%   kopt(2) mode of analysis
%        ==1 moving sum
%        ==2 moving average
%        ==3 moving median
%        ==4 frequency of events
%   kopt(3) pretreatment of data
%        ==1 none -- work with original data
%        ==2 subtract long term mean(s), i.e., work with departure series
%        ==3 convert to percent of long term mean
%        ==4 standardized departures
%   kopt(4) threshold direction (N/A unless kopt(2)==4)
%        ==1 event defined as value below threshold
%        ==2 event defined as value above threshold
%   kopt(5) threshold determination method
%        ==1 mean for period YRS(2,:);
%        ==2 median ...
%        ==3 value passed as argument thresh, which otherwise is []
%   kopt(6) how to get year vector for input time series
%        ==1 year is first col in the data matrices
%        ==2 year is stored in cv yr, read in from same .mat file as the data
%
% wsize1 (1 x 3)i  specifies window sizes to use in constructing color map 
%    col 1: start
%    col 2: increment
%    col 3: end  .  Example: [1:1:30] in annual mode means window sizes 1 to 30 yr
%           while [1:1:360] in monthly mode means window sizes 1 to 360 months
% wsize2 (1 x ?);   window sizes for which output data are saved. Example: [1 3 5 10 20 30]
%
% textin {}  
%  Series Name, to go in figure title -- 'N San Pedro'
%  Data type -- 'PPT', 'TMP'
%  Units  -- 'in', '^oF'
%  Color scheme -- 'jet'
%
% seasons{}c  names of seasons. Examples: {'NV-AP', 'MY-OC'}, {'Water Year'}
%     Not applicable for monthly analysis (kopt(1)==3) -- can set as [];
%  
% thresh(1 x ?)r  threshold for defining event in frequency analysis. [] unless
%    kopt(5)==3. ? is the number of seasons or months
% contr {} 2-element cell with plot contour control info
%  {1} (1 x ?)r contours at these values
%  {2} (1 x 2)r caxis setting for mapping color to contours
%
%*********************  FUNCTIONS CALLED
%
% movsum1.m -- moving sum
% movmed1.m -- moving median

%******************  CHECK INPUT ********************

nfiles=size(filesin,2); % Number of input files

% --------------- kopt
if length(kopt)~=6;
   error('kopt should be 1 x 6');
end
if  kopt(1)<1 | kopt(1)>3;
   error('kopt(1) must be 1, 2, or 3');
end
if  kopt(2)<1 | kopt(1)>4;
   error('kopt(2) must be 1, 2,3 or 4');
end
if  kopt(3)<1 | kopt(3)>4;
   error('kopt(3) must be 1, 2,3, or 4');
end
if  kopt(4)<1 | kopt(4)>2;
   error('kopt(4) must be 1 or 2');
end
if  kopt(5)<1 | kopt(5)>3;
   error('kopt(5) must be 1, 2, or 3');
end
if kopt(6)<1 | kopt(6)>2;
   error('kopt(6) must be 1 or 2');
end


if kopt(1)==3;  % Monthly mode
   if length(jcol)~=1;
      error('Monthly mode, jcol must be scalar')
   end
   if nfiles~=1;
      error('Monthly mode requires that all data from one file (nfiles must be 1)');
   end
   if ~isempty(seasons);
      error('Monthly mode should have seasons empty');
   end
	if kopt(3)==1;
		warndlg('Moving whatever over over monthly data is month dependent');
	end
   nseas=12; % number of 'seasons'
   timeunit='Month'; % time unit used in labeling
elseif kopt(1)==2; % multi season mode
   if length(jcol)<2;
      error('multi-season mode, but length of jcol less than 2');
   end
   nseas=length(jcol); % number of seasons
   if size(seasons,2)~=nseas;
      error('jcol says 2 seasons, but cell array seasons does not have 2 season names');
   end
   if ~(nfiles==1 | nfiles==nseas);
      error('Multi-season mode requires that nfiles be 1 or equal to number of seasons');
   end
   
else ; % kopt(1) is 1 -- annual mode
   if nfiles~=1; 
      error('Annual mode, but more than one file in filesin');
   end
   if length(jcol)~=1;
      error('Annual mode requires that jcol be length 1');
   end
   if size(seasons,2)~=1;
      error('Annual mode requires that cell array seasons name just one season');
   end
   nseas=1;
end; % if kopt(1)==3


%----------------  YRS

[mtemp,ntemp]=size(YRS);
if ~(mtemp==3 & ntemp==2);
   error('YRS should be 3 x 2');
end

%--------------    Windows -- wsize1 and wsize2
[mtemp,ntemp]=size(wsize1);
if mtemp~=1 | ntemp~=3;
   error('wsize1 should be 1 x 3');
end
[mtemp,ntemp]=size(wsize2);
if mtemp~=1 | ntemp<2;
   error('wsize2 should be rv of length at least 2');
end
nsave=length(wsize2);    % number of slabs for which info will be saved
% Build windows
wind1=[wsize1(1):wsize1(2):wsize1(3)]; % computation will be done for these
		% window sizes
wind2=wsize2; % Selected detailed info saved for these window sizes

%---------------  thresh
if kopt(5)==3; % thrhold to be passed as input argument, not empty
   if isempty(thresh);
      error('kopt(5) is 3, but thresh is empty');
   else;
      [mtemp,ntemp]=size(thresh);
      if mtemp~=1 | ntemp~=nseas;
         error('thresh must be a rv of size nseas')
      end
   end
else
   if ~isempty(thresh);
      error('kopt(5) is 3, but thresh is not empty');
   end
end


%----------------  textin
ntemp=size(textin,2);
if ntemp~=4;
   error('textin should be cell array of length 4');
end
nameser=char(textin(1));  % series name, to go on plot
if ~any(strcmp(textin(2),{'PPT','TMP'}));
   error('Invalid entry in textin for data type');
else
   dtype=char(textin(2));
end
if ~any(strcmp(textin(3),{'in','oF','mm'}));
   error('Invalid entry of data units for textin(3)');
else
   units=char(textin(3));
end
if ~any(strcmp(textin(4),{'jet','cool'}));
   error('Invalid entry for colorscheme in textin(4)');
else
   colorscm=char(textin(4));
end


%-------------------  contr
ntemp=size(contr,2);
if ntemp~=2;
   error('contr must be cell array of length 2');
end
v=contr{1};  % specifies contours, or number of contours if a scalar
clim1=contr{2}; % these data values to be mapped to extremes of the color map



%*****************   BUILD THE TIME SERIES

% Initialize
a=NaN;
yr1=(YRS(1,1):YRS(1,2))';  % year vector for plot period
yr2=(YRS(2,1):YRS(2,2))';  % year vector for period to be used for 'standardizing' stats, 
yr3=(YRS(3,1):YRS(3,2))';  % year vector for complete analysis period
	% such as the 'long-term' mean
nyrs1=length(yr1); % length of plot period
nyrs2=length(yr2); % length of means period
nyrs3=length(yr3); % length of computation period


%-----------  Compute year vectors, offset indices, and logical vectors
%  needed to deal with padding of beginnning of data and with getting
%  the correct data back after running analysis on the padded data.
%  Data will be padded with the long-term means.  The number of years to
%  pad will be set equal to the equivalent number of years in the widest
%  window specified by wind1
%
% If annual or single-season mode, 
% series, will pad with maximum window width of years.  If multi-season,
% assume that nseas seasons in a year and pad with max(wind1)/nseas.
% If monthly, know that wind1 entries are numbers of months.  If
% max(wind1) is not a multiple of 12, will pad with number of 
% years equal to the next highest multiple of 12 months.
%
%

% Compute number of years to pad, and offset indices
maxwind=max(wind1);
if kopt(1)==1; % Annual or single-season mode
	nyrpad=maxwind; % pad with this many years
	i1=nyrpad+1;  % index to start year of computation period
	i2=i1+yr1(1)-yr3(1); % index to start year  of plot period
	i1a=i1; % index in terms of observations rather than years
   i2a=i2; % ditto
   isp1=nyrs3+nyrpad;
elseif kopt(1)==2; % Multi-season mode
	if rem(maxwind,nseas)~=0;
		error('maximum of wind1 must be multiple of nseas');
	else
		nyrpad=maxwind/nseas;
		i1=nyrpad+1;
		i2=i1+yr1(1)-yr3(1);
      i1a=nyrpad*nseas+1;
      i2a=i1a + nseas*(yr1(1)-yr3(1));
      isp1=(nyrs3+nyrpad)*nseas;
	end
else; % kopt(1)==3 : monthly mode
	temp1=rem(maxwind,12)~=0;
	if temp1~=0;
		nyrpad=ceil(maxwind/12);
	else
		nyrpad=maxwind/12;
	end
	i1=nyrpad+1;
	i2=i1+yr1(1)-yr3(1);
	i1a=nyrpad*nseas+1;
   i2a=i1a + nseas*(yr1(1)-yr3(1));
   isp1=(nyrs3+nyrpad)*nseas;
end	



%---------- Annual Mode

if kopt(1)==1;
   clear X Y Z
   pf1=char(filesin(1)); % path and filename of file holding the time series
   % Get the .mat file that has matrix
   j=jcol-1; % the col of the matrix in pf1 that will have the desired data,
   % after year col is stripped off
   eval(['load ' pf1]);
   % Figure out if data file is X or Y, or Z, and make it X
   if exist('X')~=1;
      if exist('Y')==1;
         X=Y;
      elseif exist('Z')==1;
         X=Z;
      else
         error('Input .mat file has no matrix by name of X,Y or Z');
      end
   end
   
   %Get year vector for X
   [mX,nX]=size(X);
   yr = X(:,1);
   if nX<jcol;
      error('Too few cols in X for a year and specified jcol entry');
   end
   
   % Lop year col off X
   X(:,1)=[];
   
   % Check out YRS
   yrgo1=yr(1);
   yrsp1=yr(mX);
   if any(YRS(:,1)<yrgo1);
      clc
      disp(char(pf1));
      error('Starting years in YRS inconsistent with years in data matrix');
   end
   if any(YRS(:,2)>yrsp1);
      clc
      disp(char(pf1));
      error('End years in YRS inconsistent with years in data matrix');
   end
   
   % Pull data for analysis period 
   L1= yr>=YRS(3,1) & yr<=YRS(3,2); % years to plot
   z1=X(L1,j);
   
   % Pull data for stats period 
   L2= yr>=YRS(2,1) & yr<=YRS(2,2); % years to plot
   z2=X(L2,j);
   
   % Check that no NaN data in analysis period 
   if any(isnan(z1));
      error('NaN data in period for analysis and plotting');
   end

   % Get monthly 'long term' mean and std dev
   xmean=nanmean(z2);
   xstd = nanstd(z2);
	 xmed=median(z2);
   
   nobs=length(z1); % number of "observations" in the analysis period
   tobs=(1:nobs)'; % sequential 'time' vector for observations in analysis pd
	 nobsplt=nseas*nyrs1;
	 tobsplt=(1:nobsplt)';


	 %  Pad front of data with long-term means
	 xpad=repmat(xmean,nyrpad,1);
	 z1=[xpad; z1];


   % Pretreat data
   if kopt(3)==1; % use original data
      z3=z1;
   elseif kopt(3)==2; % departures from mean
      z3=z1-xmean;
   elseif kopt(3)==3; % pctg of normal
      z3=100*z1/xmean;
   else; % standardized departures
      z3=(z1-xmean)/xstd;
   end

	 % Convert z3 to event series, if desired
	 if kopt(2)==4;
	 	% Get or compute threshold
     if kopt(5)==1; % use mean from stats period
            zcrit=xmean;
     elseif kopt(5)==2; % use median
            zcrit=xmed;
     elseif kopt(5)==3; % use value in thresh
            zcrit=thresh;
     else; % reserved
     end
     % Mark event, depending on whether below or above critical value
     if kopt(4)==1; % below is event
     	z3=z3<zcrit;
     else
     	z3=z3>zcrit;
		end
	 end; % if kopt(2)==4

elseif kopt(1)==3; % monthly mode

   clear X Y Z
   pf1=char(filesin(1)); % path and filename of file holding the time series
   % Get the .mat file that has matrix
   j=(2:13) + (jcol-1)*12; % the cols holding desired Jan-Dec data,
   eval(['load ' pf1]);
   % Figure out if data file is X or Y, or Z.  Put the data in
	% X, with year in col 1
   if exist('X')~=1;
      if exist('Y')==1;
         X=[Y(:,1) Y(:,j)];
      elseif exist('Z')==1;
         X=[Z(:,1) Z(:,j)];
      else
         error('Input .mat file has no matrix by name of X,Y or Z');
      end
	else; % X exists
		X=[X(:,1) X(:,j)];
   end
   
   %Get year vector for X
   [mX,nX]=size(X);
   yr = X(:,1);
   if nX~=13;
      error('X should have 13 cols at this point');
   end
   
   % Lop year col off X
   X(:,1)=[];
   
   % Check out YRS
   yrgo1=yr(1);
   yrsp1=yr(mX);
   if any(YRS(:,1)<yrgo1);
      clc
      disp(char(pf1));
      error('Starting years in YRS inconsistent with years in data matrix');
   end
   if any(YRS(:,2)>yrsp1);
      clc
      disp(char(pf1));
      error('End years in YRS inconsistent with years in data matrix');
   end
   
   % Pull block of monthly data for analysis period 
   L1= yr>=YRS(3,1) & yr<=YRS(3,2); % years to analyze
   Z1=X(L1,:);
   
   % Pull data for stats period 
   L2= yr>=YRS(2,1) & yr<=YRS(2,2); % years for mean, std dev
   Z2=X(L2,:);
   
   % Check that no NaN data in analysis period 
   if any(any(isnan(Z1)));
      error('NaN data in period for analysis');
   end

   % Get rv's of 12 monthly 'long term' means, medians and std devs
   xmean=nanmean(Z2);
   xstd = nanstd(Z2);
	 xmed = nanmedian(Z2);


	 % Pad front of data with long-term monthly means
	 xpad=repmat(xmean,nyrpad,1);
	 Z1=[xpad; Z1];
   mZ1=size(Z1,1); % number of years, including padded section

	 % Dupe rows to make matrices of means and std devs, same row size as Z1 
	 Xmean=repmat(xmean,mZ1,1);
	 Xstd=repmat(xstd,mZ1,1);
	 Xmed=repmat(xmed,mZ1,1);
   
   nobs=12*nyrs3; % number of "observations" in the analysis period
   tobs=(1:nobs)'; % sequential 'time' vector for observations
	 nobsplt=12*nyrs1;
	 tobsplt=(1:nobsplt)';

   % Pretreat data
   if kopt(3)==1; % use original data
      Z3=Z1;
   elseif kopt(3)==2; % departures from mean
      Z3=Z1 - Xmean;
   elseif kopt(3)==3; % pctg of normal
      Z3=100*(Z1 ./ Xmean);
   else; % standardized departures
      Z3=(Z1-Xmean) ./ Xstd;
   end


	 % If Desired, convert Z3 into an 'event' series - 0 or 1
	 if kopt(2)==4;
 	 	% Get or compute thresholds for the nseas seasons or 12 months -- a rv
	 	% with nseas elements
     if kopt(5)==1; % use means from stats period
     	zcrit=Xmean;
     elseif kopt(5)==2; % use median
        zcrit=Xmed; 
     elseif kopt(5)==3; % use value in thresh
        zcrit=repmat(thresh,mZ1,1);
		else; % reserved
     end
		
		% Mark event, depending on whether below or above critical value
     if kopt(4)==1; % below is event
     	Z3=(Z3<zcrit);
     else
			Z3=(Z3>zcrit);
     end
	 end; % if kopt(2)==4;
		
   % String monthly matrix into a single time series
	 Z3=Z3';
   z3=Z3(:);
   
   %-------------------------------------------------
   
   
elseif kopt(1)==2; % multi season mode
	disp('Looping over seasonal files to build data matrices');
	for n =1:nfiles;
		clear X Y Z
		pf1=cellstr(filesin(1,n)); % path and filename
		seasnm=cellstr(seasons(1,n)); % name of season -- used in labeling
		% Get the .mat file that has matrix
		 j=jcol(n)-1; % the col of the matrix in pf1 that will have the desired data,
   			% after neglecting the year col
		disp([  '   Season ' char(seasnm) ' from ' char(pf1) ', variable # ' int2str(j)]);
		eval(['load ' char(pf1)]);

		% Figure out if data file is X or Y, or Z, and make it X
		if exist('X')~=1;
			if exist('Y')==1;
				X=Y;
			elseif exist('Z')==1;
				X=Z;
			else
			error('Input .mat file has no matrix by name of X,Y or Z');
			end
		end

		% Get year vector for X
		[mX,nX]=size(X);
		yr = X(:,1);
		if nX<jcol(n);
			error('Too few cols in X for a year and specified jcols entry');
		end
   
		% Lop year col off X
		X(:,1)=[];

		% Check out YRS
		yrgo1=yr(1);
		yrsp1=yr(mX);
		if any(YRS(:,1)<yrgo1);
			clc
			disp(char(pf1));
			error('Starting years in YRS inconsistent with years in data matrix');
		end
		if any(YRS(:,2)>yrsp1);
			clc
			disp(char(pf1));
			error('End years in YRS inconsistent with years in data matrix');
		end
   
		% Pull data for analysis period and put in col of Z1
		L1= yr>=YRS(3,1) & yr<=YRS(3,2); % years to plot
		Z1(:,n)=X(L1,j);

		% Pull data for stats period and put in col of Z2
		L2= yr>=YRS(2,1) & yr<=YRS(2,2); % years to plot
		Z2(:,n)=X(L2,j);

	end; % for n=1:nfiles
	
	% Check that no NaN data in analysis period for any series
	if any(any(isnan(Z1)));
  	error('NaN data in period for analysis');
	end

	% Get rv's of nseas seasonal 'long term' means, medians and std devs
  xmean=nanmean(Z2);
  xstd = nanstd(Z2);
	xmed = nanmedian(Z2);

  % Pad front of data with long-term monthly means
	xpad=repmat(xmean,nyrpad,1);
	Z1=[xpad; Z1];
  mZ1=size(Z1,1); % number of years, including padded section

	% Dupe rows to make matrices of means and std devs, same row size as Z1 
	 Xmean=repmat(xmean,mZ1,1);
	 Xstd=repmat(xstd,mZ1,1);
	 Xmed=repmat(xmed,mZ1,1);

   nobs=nseas*nyrs3; % number of "observations" in the analysis period
   tobs=(1:nobs)'; % sequential 'time' vector for observations
	 nobsplt=nseas*nyrs1;
	 tobsplt=(1:nobsplt)';


  % Pretreat data
  if kopt(3)==1; % use original data
     Z3=Z1;
  elseif kopt(3)==2; % departures from mean
     Z3=Z1 - Xmean;
  elseif kopt(3)==3; % pctg of normal
     Z3=100*(Z1 ./ Xmean);
  else; % standardized departures
     Z3=(Z1-Xmean) ./ Xstd;
  end


	 % If Desired, convert Z3 into an 'event' series - 0 or 1
	 if kopt(2)==4;
 	 	% Get or compute thresholds for the nseas seasons or 12 months -- a rv
	 	% with nseas elements
     if kopt(5)==1; % use means from stats period
     	zcrit=Xmean;
     elseif kopt(5)==2; % use median
        zcrit=Xmed; 
     elseif kopt(5)==3; % use value in thresh
        zcrit=repmat(thresh,mZ1,1);
		else; % reserved
     end
		
		% Mark event, depending on whether below or above critical value
     if kopt(4)==1; % below is event
     	Z3=(Z3<zcrit);
     else
			Z3=(Z3>zcrit);
     end
	 end; % if kopt(2)==4;
		
   % String multi-season matrix into a single time series
	 Z3=Z3';
   z3=Z3(:);

end; % if kopt(1)==3;



%***************** ANALYZE THE SERIES -- COMPUTE THE MATRIX FOR THE COLOR PLOT

% Initialize data matrices
colmax=length(wind1); 
D=a(ones(nyrs3*nseas,1),ones(colmax,1)); % will hold full computation period

% Initialize the color map data matrix. Note that for an n-yr moving ave, need
% n-1 leading values to get a filtered data to plot at time 1.  I will assume the
% leading values are equal to the long-term mean.  For monthly data, must
% retain individual monthly means if original data used. Likewise for 
% multi-season data

if kopt(1)==1; % annual mode
   kk1=[1 1]; % options for movsum1: moving period plotted at end, and omit plot 
   % loop over moving sums
   for n =1:length(wind1);
      nn=wind1(n); % current window width (in yrs)
		 igo1 =i1a-nn+1; % start index of analysis data
		 z4=z3(igo1:isp1);  % raw data, before smoothing
		 tt = (1: (isp1-igo1+1))'; %  relative time index for z4
      
      if kopt(2)==1 | kopt(2)==2 | kopt(2)==4; % moving sum, ave, or freq
         [y4,tt4]=movsum1(z4,tt,nn,kk1);
         if kopt(2)==2;
            y4=y4/nn; % mov sum to mov ave
         end
      elseif kopt(2)==3; % moving median
         [y4,tt4]=movmed1(z4,tt,nn,1);
      end

            
      % Check length of y4
      if length(y4)~=nyrs3;
         error('Returned moving ave vector not of correct length');
      end
         
      
      % Put moving whatever series in col of D
      D(:,n)=y4;
   end; % for n=1:length(wind1)

elseif kopt(1)==3 | kopt(1)==2;; % monthly mode or multi season mode
  kk1=[1 1]; % options for movsum1: mov period plotted at end, omit plot 
   
	% loop over moving sums
   for n =1:length(wind1);
      nn=wind1(n);
		 igo1 =i1a-nn+1; % start index of analysis data
		 z4=z3(igo1:isp1);  % raw data, before smoothing
		 tt = (1: (isp1-igo1+1))'; %  relative time index for z4
      if kopt(2)==1 | kopt(2)==2 | kopt(2)==4; % moving sum or moving ave or
				% freq of events
         [y4,tt4]=movsum1(z4,tt,nn,kk1);
         if kopt(2)==2;
            y4=y4/nn; % mov sum to mov ave
         end
      elseif kopt(2)==3; % moving median
         [y4,tt4]=movmed1(z4,tt,nn,1);

      end

      % Check length of y4
      if length(y4)~=nobs;
         error('Returned moving ave vector not of correct length');
      end
         
      
      % Put moving whatever series in col of D
      D(:,n)=y4;
   end; % for n=1:length(wind1)
end; % if kopt(1)==3



%******************** CULL PORTION OF D TO BE USED FOR THIS PLOT

% D now has data for full YRS(3,:) period.  Want only that for YRS(1,:);


if kopt(1)==1; % annual or single season mode
	L1temp=yr3>=YRS(1,1) & yr3<=YRS(1,2);
	D=D(L1temp,:);
else;  % multiseason or monthly mode
	istart=i2a-i1a+1;
	iend=istart +   nyrs1*nseas -1;
	D=D(istart:iend,:);
end



%****************** CULL OUTPUT DATA TO BE RETURNED -- subset of rows of D

% Expand rv wind2 to matrix, same row size as wind1
W2=repmat(wind2,colmax,1);

% Make comparison matrix of wind1
W1=wind1'; % wind1 to cv
W1=repmat(W1,1,length(wind2));

D1=D(any((W1==W2)'),:);



%******************  PLOT THE FIGURE, ALLOWING FOR INTERACTIVE CONTOUR LABELING, ETC

% Note that D now has years (or months or seasons) going down col, 
% length of window length to right

% in plotting D, I expect time on x axis, and longest smoothing period
% at top of figure.  Rotate matrix to proper position
D=rot90(D,1);  % ccw rotation, 90 degrees
D=flipud(D);
D1=rot90(D,1);
D1=flipud(D1);

% Get plotting variable for time
if kopt(1)==1; % annual mode
	tplot=yr1;
else;
	tplot=yr1(1) + (tobsplt-1)/nseas;
end

% Make color map
figure(1)
[C,H]=contourf(tplot,wind1,D,v);
if isempty(clim1),
   clim1=[min(min(D)) max(max(D))];
end

caxis(clim1);
set(gcf,'ColorMap',flipud(jet));
%clabel(C,H);

%------------  Build variable labeling info for title and axes

tit1=nameser;
tit3=dtype;
tit4=sprintf(', %4.0f-%4.0f',YRS(1,:));
tit10=[' Contours (' units '): '];
tit11=sprintf('%4.1f ',v);

if kopt(1)==1; % single season mode
   tit2='Seasonal ';
   tit7=[' / Season: ' char(seasons(1))]; 
   ystr1='yr';
elseif kopt(1)==2; % multi season mode
   tit2='Seasonal ';
   % Make a single string out of the seasons
   tit7=' / Seasons: ';
   for m=1:size(seasons,2);
      str1=[char(seasons(m)) ' '];
      tit7=[tit7 str1];
   end
   ystr1='seasons';
else ; % monthly mode
   tit2='Monthly ';
   tit7=' ';
   ystr1='months';
end; % if kopt(1)==1



if kopt(2)==1;
   tit8='Moving Sum; ';
   tit9='';
elseif kopt(2)==2;
   tit8='Moving Average; ';
   tit9='';
elseif kopt(2)==3;
   tit8='Moving Median; ';
   tit9='';
else; % kopt(2)==4;
   tit8='Moving Frequency of Events: x';
   if kopt(4)==1; % below thresh
      tit9a='<';
          
   else
      tit9a='>';
      set(gcf,'ColorMap',jet);
   end
   tit9b=num2str(zcrit(1,:));
   tit9=[tit9a tit9b ';'];
end; % if kopt(2)==1


if kopt(3)==1; % original data
   tit5='Original';
   tit6='';
elseif kopt(3)==2;
   tit6=sprintf('%4.0f-%4.0f Mean',YRS(2,:));
   tit5='Departures from ';
elseif  kopt(3)==3;
   tit6=sprintf('%4.0f-%4.0f Mean',YRS(2,:));
   tit5='Percentage of ';
else; % stdzd anomalies
   tit6=sprintf('%4.0f-%4.0f Mean',YRS(2,:));
   tit5='Stdzd anomalies from ';
end


tline1=[tit1 '; Moving Window Summary of ' tit2 tit3 tit4];
tline2=['Data: ' tit5 tit6 tit7];
tline3=['Method: ' tit8 tit9 tit10 tit11];


% Change position to allow room for 3-line title
set(gca,'Position',[.13 .11 .775 .75]);

% Add title and axes labels
title({tline1, tline2, tline3});
xlabel('Ending Year of Windowed Period');
ylabel(['Window Size (' ystr1 ')']);

% Add color bar
colorbar

% Make ticks point out
set(gca,'TickDir','out');

  % Add grid lines for x axis
      set(gca,'Layer','top','XGrid','on');


% Manually label contours
clabel(C,H,'manual')

E=D(:);
% Interactively decide on contours 
k1=1;
answer3=clim1;
answer2=v;
while k1==1;
   answer1=questdlg('Happy with the contours and color?');
   switch answer1
   case 'No'
      iprompt={'Enter desired contour values'}
      ititle='CHANGING THE CONTOURS';
      %idef=prctile(E,[20 40 60 80]); 
      idef2=answer2;
      ilineno=1;
      answer2=dlgi001(ititle,iprompt,ilineno,idef2);
      
      iprompt={'Enter desired caxis values'}
      ititle='CHANGING THE CAXIS';
      idef3=answer3;
      ilineno=1;
      answer3=dlgi001(ititle,iprompt,ilineno,idef3);

      [C,H]=contourf(tplot,wind1,D,answer2);
      caxis(answer3);
      %clabel(C,H,answer2);
      
      % Change position to allow room for 3-line title
      set(gca,'Position',[.13 .11 .775 .75]);

      % Add title and axes labels
      tit11=sprintf('%4.1f ',v);
      tline3=['Method: ' tit8 tit9 tit10 tit11];
      title({tline1, tline2, tline3});
      xlabel('Ending Year of Windowed Period');
      ylabel(['Window Size (' ystr1 ')']);

      % Add color bar
      colorbar

      % Make ticks point out
      set(gca,'TickDir','out');
      
      % Add grid lines for x axis
      set(gca,'Layer','top','XGrid','on');
      

      % Manually label contours
      clabel(C,H,'manual')

      
      
   case 'Cancel'
      fclose('all');
      return
   case 'Yes'
      k1=0;
   end
end



%*********************  SUMMARY TEXT FILE

[mD,nD]=size(D);

% Rank D along rows
D=D';
[S,I]=sort(D);  % sorted ascending

% If frequency of events mode, will be interested in highest values
if kopt(2)==4;
   S=flipud(S); % sorted descending
   I=flipud(I);  % ditto
end


% Compute year, month corresponding to index I
% Convert sequential month of end of runs to year,month
Iend=I;  % ending sequential obs of period
monend=rem(Iend,nseas);
L1temp=monend==0;
yearend=yr1(1)+floor(Iend/nseas);
monend(L1temp)=nseas;
yearend(L1temp)=yr1(1)+Iend(L1temp)/nseas-1;

[file2,path2]=uiputfile('*.txt','File for summary of moving window analysis');
pf2=[path2 file2];
fid2=fopen(pf2,'w');

%-----------------  Opening text
fprintf(fid2,'%s\n\n','MOVING WINDOW SUMMARY');
fprintf(fid2,'%s\n',tline1);
fprintf(fid2,'%s\n',tline2);
fprintf(fid2,'%s\n\n',tline3);


fprintf(fid2,'Units of window width = %s\n',ystr1);
fprintf(fid2,'Units of values = %s\n',units);



%---- Table of Lowest value for each window size



% Loop over windows, from narrowest to widest

fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));

fprintf(fid2,'%s\n\n','SUMMARY OF LOWEST VALUE FOR EACH WINDOW SIZE');

head2a='Window    End                 ';
head2b='Width     Year          Value  ';
fprintf(fid2,'%s\n',head2a);
fprintf(fid2,'%s\n\n',head2b);

for n =1:mD;
   howwide=wind1(n);
   str11=sprintf('%2.0f',n);
   str12=sprintf('%3.0f\t     ',howwide);
   str13=sprintf('%2.0f/',monend(1,n));
   str14=sprintf('%4.0f\t',yearend(1,n));
   str15=sprintf('%g',S(1,n));
   
   if nseas>1;
      strall1=[ str12  str13 str14  str15];
   else
      strall1=[ str12  str14   str15];
   end
   
   fprintf(fid2,'%s\n',strall1);
   
end


%------------  Table of highest values

fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));

fprintf(fid2,'%s\n\n','SUMMARY OF HIGHEST VALUE FOR EACH WINDOW SIZE')
head2a='Window    End                 ';
head2b='Width     Year          Value  ';
fprintf(fid2,'%s\n',head2a);
fprintf(fid2,'%s\n\n',head2b);

for n =1:mD;
   howwide=wind1(n);
   str11=sprintf('%2.0f',n);
   str12=sprintf('%3.0f\t     ',howwide);
   str13=sprintf('%2.0f/',monend(nD,n));
   str14=sprintf('%4.0f\t',yearend(nD,n));
   str15=sprintf('%g',S(nD,n));
   
   if nseas>1;
      strall1=[ str12  str13 str14  str15];
   else
      strall1=[ str12  str14   str15];
   end
   
   fprintf(fid2,'%s\n',strall1);
   
end


% Tables of selected ranked events

nranked=10;  % hard coded for now -- will make tables of the 10-worst moving events
% in each window specified by wind2.
% save printin S wind1 wind2 nranked fid2 yearend monend
%print01(S,wind1,wind2,nranked,fid2,yearend,monend);




%--------------- Cull desired window sizes

% Expand rv wind2 to matrix, same row size as wind1
W2=repmat(wind2,size(wind1,2),1);

% Make comparison matrix of wind1
W1=wind1'; % wind1 to cv
W1=repmat(W1,1,length(wind2));


% Logical pointer to desired windows
Ltemp=any((W1==W2)');  % a rv
wsizes=wind1(Ltemp); % desired window sizes

% Cull desired info on moving whatever, end year and end season/month
S1=S(1:nranked,Ltemp);
yearend=yearend(1:nranked,Ltemp);
monend=monend(1:nranked,Ltemp);

% Separate from previous printout
fprintf(fid2,'%s\n\n\n\n',blanks(5));

fmt1='%4.0f  %4.0f/%2.0f\t\t %g\n';
head2='Rank  End(Yr/mo)          Value';

%-------------  LOOP OVER SELECTED WINDOW SIZES
for n1 = 1:length(wsizes);
   width1=wsizes(n1);  % window width
   head1=sprintf('\n\n\n WINDOW SIZE  = %4.0f\n\n',width1);
   fprintf(fid2,'%s',head1);
   fprintf(fid2,'%s\n\n',head2);
   
   % Loop over ranks
   for n2=1:nranked;
      fprintf(fid2,fmt1,n2,yearend(n2,n1),monend(n2,n1),S1(n2,n1));
   end
  
      
end



fclose(fid2);



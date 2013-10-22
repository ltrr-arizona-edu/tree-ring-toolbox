function D=moving1(filesin,jcol,YRS,kopt,wsize1,wsize2,textin,seasons,thresh,contr);
% moving1:  color map of moving sum, average, median, or frequency of events
% CALL: D=moving1(filesin,jcol,YRS,kopt,wsize1,wsize2,textin,seasons,thresh,contr);
%
% Meko 4-30-97
%
%******************* IN *****************************************************
%
% filesin {}  cell array, 1 row, containing paths/names of files holding input
%   time series data.  Col size of filesin tells how many files hold input data
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
%   row 1: period for analysis and plot
%   row 2: period to compute means on (anomalies might be needed)
%
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
% thresh(1 x 1)r  threshold for defining event in frequency analysis. [] unless
%    kopt(5)==3
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
if ~(mtemp==2 & ntemp==2);
   error('YRS should be 2 x 2');
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


%---------------  thresh
if kopt(5)==3; % thrhold to be passed as input argument, not empty
   if isempty(thresh);
      error('kopt(5) is 3, but thresh is empty');
   else;
      [mtemp,ntemp]=size(thresh);
      if mtemp~=1 | ntemp~=1;
         error('thresh must be a scalar')
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
if ~any(strcmp(textin(3),{'in','oF'}));
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
v=contr{1};  
clim1=contr{2};



%*****************   BUILD THE TIME SERIES

% Initialize
a=NaN;
yr1=(YRS(1,1):YRS(1,2))';  % year vector for plot period, which is also the analysis period
yr2=(YRS(2,1):YRS(2,2))';  % year vector for period to be used for 'standardizing' stats, 
% such as the 'long-term' mean
nyrs1=length(yr1); % length of analysis & plot period
nyrs2=length(yr2); % length of means period

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
   
   % Pull data for plot period 
   L1= yr>=YRS(1,1) & yr<=YRS(1,2); % years to plot
   z1=X(L1,j);
   
   % Pull data for stats period 
   L2= yr>=YRS(2,1) & yr<=YRS(2,2); % years to plot
   z2=X(L2,j);
   
   % Check that no NaN data in plot period 
   if any(isnan(z1));
      error('NaN data in period for analysis and plotting');
   end

   % Get monthly 'long term' mean and std dev
   xmean=nanmean(z2);
   xstd = nanstd(z2);
   
   nobs=length(z1); % number of "observations" in the plot period
   tobs=(1:nobs)'; % sequential 'time' vector for observations

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
   
   % Pull block of monthly data for plot period 
   L1= yr>=YRS(1,1) & yr<=YRS(1,2); % years to plot
   Z1=X(L1,:);
   
   % Pull data for stats period 
   L2= yr>=YRS(2,1) & yr<=YRS(2,2); % years to plot
   Z2=X(L2,:);
   
   % Check that no NaN data in plot period 
   if any(any(isnan(Z1)));
      error('NaN data in period for analysis and plotting');
   end

   % Get rv's of 12 monthly 'long term' means, medians and std devs
   xmean=nanmean(Z2);
   xstd = nanstd(Z2);
	xmed = nanmedian(Z2);

	% Dupe rows to make matrices of means and std devs 
	Xmean=repmat(xmean,nyrs1,1);
	Xstd=repmat(xstd,nyrs1,1);
   
   nobs=12*size(Z1,1); % number of "observations" in the plot period
   tobs=(1:nobs)'; % sequential 'time' vector for observations

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
   
		% Pull data for plot period and put in col of Z1
		L1= yr>=YRS(1,1) & yr<=YRS(1,2); % years to plot
		Z1(:,n)=X(L1,j);

		% Pull data for stats period and put in col of Z2
		L2= yr>=YRS(2,1) & yr<=YRS(2,2); % years to plot
		Z2(:,n)=X(L2,j);

	end; % for n=1:nfiles
	
	% Check that no NaN data in plot period for any series
	if any(any(isnan(Z1)));
  	error('NaN data in period for analysis and plotting');
	end

	% Get rv's of nseas seasonal 'long term' means, medians and std devs
  xmean=nanmean(Z2);
  xstd = nanstd(Z2);
	xmed = nanmedian(Z2);

	% Dupe rows to make matrices of means and std devs 
	Xmean=repmat(xmean,nyrs1,1);
	Xstd=repmat(xstd,nyrs1,1);

	% Compute number of "observations", which is seasons times years 
	nobs=nseas*size(Z1,1); % number of "observations" in the plot period
  tobs=(1:nobs)'; % sequential 'time' vector for observations

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


	% String seasonal matrix into a single time series
	Z3=Z3';
  z3=Z3(:);

end; % if kopt(1)==3;





%***************** ANALYZE THE SERIES -- COMPUTE THE MATRIX FOR THE COLOR PLOT

% Build windows
wind1=[wsize1(1):wsize1(2):wsize1(3)];
wind2=wsize2;

% Initialize data matrix to be contoured
colmax=length(wind1); 
D=a(ones(nyrs1*nseas,1),ones(colmax,1));

% Initialize the color map data matrix. Note that for an n-yr moving ave, need
% n-1 leading values to get a filtered data to plot at time 1.  I will assume the
% leading values are equal to the long-term mean.  For monthly data, must
% retain individual monthly means if original data used. Likewise for 
% multi-season data

if kopt(1)==1; % annual mode
   mean1=mean(z3); % will use this value as leading values
   kk1=[1 1]; % options for movsum1: moving period plotted at end, and omit plot 
   % loop over moving sums
   for n =1:length(wind1);
      nn=wind1(n);
      % slap on leading values
      if nn>1;
         nslap=nn-1;
         z4=[mean1(ones(nslap,1),:) ; z3];
         tt=(1: (nyrs1+nslap))';
      else
         nslap=[];
         z4=z3;
         tt=(1:nyrs1)';
      end
      if kopt(2)==1 | kopt(2)==2; % moving sum or moving ave
         [y4,tt4]=movsum1(z4,tt,nn,kk1);
         if kopt(2)==2;
            y4=y4/nn; % mov sum to mov ave
         end
      elseif kopt(2)==3; % moving median
         [y4,tt4]=movmed1(z4,tt,nn,1);
      else; % kopt(2) must be 4, meaning freqency of events
         % Get or compute threshold
         if kopt(5)==1; % use mean from stats period
            zcrit=mean(z2);
         elseif kopt(5)==2; % use median
            zcrit=median(z2);
         elseif kopt(5)==3; % use value in thresh
            zcrit=thresh;
         else; % reserved
         end
         % Mark event, depending on whether below or above critical value
         if kopt(4)==1; % below is event
            L1=(z4<zcrit);
         else
            L1=(z4>zcrit);
         end
         % Compute moving sum of number of events, then convert to decimal fraction
         % of events in window
         [y4,tt4]=movsum1(L1,tt,nn,kk1);
         y4=y4/nn;
      end
            
      % Check length of y4
      if length(y4)~=nyrs1;
         error('Returned moving ave vector not of correct length');
      end
         
      
      % Put moving whatever series in col of D
      D(:,n)=y4;
   end; % for n=1:length(wind1)

elseif kopt(1)==3 | kopt(1)==2;; % monthly mode or multi season mode
	% Set values to use as leading values
	if kopt(3)==1; % using original data
		mean1=xmean; % the nseas long-term means
	elseif kopt(3)==2; % using departures from long term means
		mean1=0;
	elseif kopt(3)==3; % using percentage of long-term means
		mean1=100;
	elseif kopt(3)==4; % using standardized anomalies 
		mean1=0;
	end

   kk1=[1 1]; % options for movsum1: mov period plotted at end, omit plot 
   
	% loop over moving sums
   for n =1:length(wind1);
      nn=wind1(n);
      % slap on leading values
      if nn>1;
         nslap=nn-1; % number of leading values needed
			if kopt(3)==2 | kopt(3)==3 | kopt(3)==4;
				z4=[mean1(ones(nslap,1),:) ; z3];
         	tt=(1: (nobs+nslap))';
			elseif kopt(3)==1;
				% find next multiple of 12 higher than nslap
				if rem(nslap,nseas)==0;  % nslap is multiple of nseas
					nslab=nslap/nseas;
					if kopt(2)~=4 % not in frequency of events mode
						xslab=(repmat(xmean,1,nslab))';
					else; % in frequency of events model
						if kopt(5)==1; % thresholds are the nseas long-term means
							xslab=(repmat(xmean,1,nslab))';
						elseif kopt(5)==2; % thresholds are medians
							xslab=(repmat(xmed,1,nslab))';
						else; % kopt(5) is 3, and thresholds read in
							xslab=(repmat(thresh,1,nslab))';
						end
					end
						
				else; % nslap is not multiple of nseas
					nslab=ceil(nslap/nseas);


					if kopt(2)~=4 % not in frequency of events mode
						xslab=(repmat(xmean,1,nslab))';
					else; % in frequency of events model
						if kopt(5)==1; % thresholds are the nseas long-term means
							xslab=(repmat(xmean,1,nslab))';
						elseif kopt(5)==2; % thresholds are medians
							xslab=(repmat(xmed,1,nslab))';
						else; % kopt(5) is 3, and thresholds read in
							xslab=(repmat(thresh,1,nslab))';
						end
					end
					nlopp=nseas-rem(nslap,nseas); 
					xslab(1:nlopp)=[];
				end
				z4=[xslab ; z3];
				tt=(1:(nobs+nslap))';
			end
      else
         nslap=[];
		  xslab=[];
         z4=z3;
         tt=(1:nobs)';
      end

      if kopt(2)==1 | kopt(2)==2; % moving sum or moving ave
         [y4,tt4]=movsum1(z4,tt,nn,kk1);
         if kopt(2)==2;
            y4=y4/nn; % mov sum to mov ave
         end
      elseif kopt(2)==3; % moving median
         [y4,tt4]=movmed1(z4,tt,nn,1);

      else; % kopt(2) must be 4, meaning freqency of events
         % Get or compute thresholds for the nseas seasons or 12 months -- a rv nseaselements
         if kopt(5)==1; % use means from stats period
            zcrit=xmean; % nseas-element rv
         elseif kopt(5)==2; % use median
            zcrit=xmed; % nseas-element rv
         elseif kopt(5)==3; % use value in thresh
            zcrit=thresh; % nseas-element rv
         else; % reserved
         end
         % Expand the thresholds vector by repeating for each year in plot period
         zcrit=(repmat(zcrit,1,nyrs1))';  % cv, nyrs1 length
		  % Add the leading values
		  zcrit = [xslab zcrit];
		  if length(zcrit)~=length(z4);
			error('zcrit must be same length as z4 for events mode');
		  end
		  
         % Mark event, depending on whether below or above critical value
         if kopt(4)==1; % below is event
            L1=(z4<zcrit);
         else
            L1=(z4>zcrit);
         end
         % Compute moving sum of number of events, then convert to dec fraction
         % of events in window
         [y4,tt4]=movsum1(L1,tt,nn,kk1);
         y4=y4/nn;
      end
            
      % Check length of y4
      if length(y4)~=nobs;
         error('Returned moving ave vector not of correct length');
      end
         
      
      % Put moving whatever series in col of D
      D(:,n)=y4;
   end; % for n=1:length(wind1)
end; % if kopt(1)==3



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

% Get plotting variable for time
if kopt(1)==1; % annual mode
	tplot=yr1;
else;
	tplot=yr1(1) + (tobs-1)/nseas;
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
tit10=[' Contours (' units ')'];

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
   tit9b=int2str(zcrit);
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
tline3=['Method: ' tit8 tit9 tit10];


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
%S=flipud(S); % sorted descending
%I=flipud(I);  % ditto

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


fclose(fid2);



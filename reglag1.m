function reglag1
% reglag1:   distributed-lag regression to filter and scale tree-ring indices 
% CALL: reglag1
%
% Meko 1-14-99 -- modified from earlier versions
%
%****************  IN **********
%
% No input arguments
%
% User prompted to point to input control file with following line by line:
%
% nsites -- number of tree-ring sites in input tsm
% path\filename of .mat file holding the tsm of tree-ring variables
% path to .mat files holding the water balance seasonal output
% seascode -- season code of water balance file <oc_se>
% krule (1,2 or 3);
%   1=variable entry stops when adjusted R-squared fails to increase
%   2=variable entry stops when R-squared fails to increase by at least deltrsq
%   3=variable entry stops with entry of npred th predictor
% deltarsq -- threshold required increment increase in regression R-squared
% yrgo, yrsp =start, end years of period to consider for modeling
% yrstango, yrstansp = start and end years of standardizing period optionally
%   used for scaling filtered reconstructions.
% nlag -- number of + and - lags possible in model 
% kmode (1 or 2): 
%  1=exploratory regression mode.  Each of 16 water-balance variables is predictand 
%  2=reconstruction mode
% path\filename of .mat file holding cv's  iclim and Lomit.
%   **iclim tells which which of the 16 water balance variables to reconstuct. 
%   **Lomit is logical cv telling which tree sites to omit from modeling
%   **npred is integer cv telling how many predictors to allow to enter model
% path\filename of .mat file with site names and other text information on tree-ring sites 
% path\filename of .dat file with elevations (m) of tree-ring sites
% path\filename of .dat file with longs and lats (decimal deg) of sites
% path\filename of .dat file with first, last years of chronologies
%
% SAMPLE .CTL INPUT
% 71
% c:\projs\af1\sacrflow\wrktree\tsma2
% c:\projs\af1\sacrflow\outwat\
% oc_se
% 2
% .05
% 1890 1995
% 1890 1961
% 2
% 2
% c:\projs\af1\sacrflow\wrktree\climwhch
% c:\projs\af1\sacrflow\wrktree\sitenms
% c:\projs\af1\sacrflow\wrktree\eleva.dat
% c:\projs\af1\sacrflow\wrktree\crnxya.dat
% c:\projs\af1\sacrflow\wrktree\gospa.dat
%
%******************* NOTES *************
%
%
% MORE ABOUT THE INPUT...
%
% nsites -- number of tree-ring sites.
%   This refers to the number of columns (minus year column)
%   in an input tsm of tree-ring data.  Which tree-ring series are modeled depends on
%   Lomit (see Lomit).
% path\file of .mat file holding the input tsm of tree-ring variables
%   These tree-ring variables are generally AR residuals that have been generated
%   using function arwhite2.m.  And even before running arwhite2.m, I truncated the
%   usable years of tree-ring data according to minimum acceptable sample size
%   (ASS year).  Best to check out /projs/af1/sacrflow/steps.doc for sequence of
%   operations in a full-blown analysis.  The tsm must also have AR prewhitening
%   information.
% path to .mat files holding the water balance seasonal output
%   To prepare the climate data, I did lots of things. P and T had to be interpolated.
%   I ran functions that generated PDSI, Z-index and other quantities for site-centered
%   climate data.  See /projs/af1/sacrflow/steps.doc.
% seascode -- season code of water balance file <oc_se>
%   This string is used in labeling.  'oc-sep' indicates water-year total or average.
% krule (1,2 or 3);
%   1=variable entry stops when adjusted R-squared fails to increase
%   2=variable entry stops when R-squared fails to increase by at least deltrsq
%   3=variable entry stops after npred th predictor has entered (see below)
% deltarsq -- threshold required increment increase in regression R-squared. Has an
%   effect only if krule==2. Otherwise, ignored.
% yrgo, yrsp =start, end years of period to consider for modeling
%   Time coverage by climate and tree-ring data might vary from site to site. yrgo and
%   yrsp are used to bracket the ***longest*** possible calibration period that might
%   be considered for any site.  A site might not have enough valid data to calibrate
%   over that full period. Then some subset of years is used.
% yrstango, yrstansp = "common period", start and end years for scaling filtered reconstructons
%    If [NaN NaN], no scaling is desired. If scaling desired, the full-period 
%    regression model (with same predictors as selected for full period) is
%    re-fit to the common period, and the R-squared values are used to 
%    scale the reconstructions relative to one another.  Mean and std dev are computed
%    for the common period, then used to convert recon to z-scores.  Then R-sqared
%    used to scale the z-scores.  All series with Lomit~=1
%    must have sufficient data coverage to allow modeling over the common period. 
%    Scaling ensures that standard deviations of the scaled reconstruction
%    computed for the common period are equal to the R-squared values of the full-period
%    model that is re-fit to the common period.  If no scaling desired, simply set
%    yrstango,yrstansp to [ Nan NaN].
% nlag -- number of + and - lags possible in model
%    This many lags will be considered as potential predictors in fitting the models.  
%    Which lags actually enter the models depends on the data.  Variables enter only as
%    long as adjusted R-squared increases.  The settings for nlag
%    also define "n" in the "leave-n-out" cross-validation (see below)
% kmode (1 or 2): 
%  1=exploratory regression mode.  
%    Focus is deciding whether one type of water-balance variable might be preferable
%    to another as a predictand for the reconstruction models.
%    Each of 16 water-balance variables is as tried as predictand in separate models. 
%    No long-term reconstructions are generated.
%    No model validation (split-sample or cross-validation) is done
%    Only the part of the tree-ring tsm for [yrgo yrsp] is used, as you do
%    not need long-term filtered, scaled series.  
%  2=reconstruction mode
%    Long-term single-site reconstructions (filtered, scaled tree-ring indices) are
%    generated using a specified water-balance variable for each tree-ring site.  In the
%    regression, verification by both split-sample and cross-validation are done
% path\file of .mat file holding cv's  iclim and Lomit.
%   **iclim is a scalar that tells which which of the 16 water balance variables to reconstuct. This makes
%     sense only for kmode==2.  If kmode==1, all 16 w-b variables are tried as predictand.
%     iclim is an integer cross-referenced to type of variable:
%      1 -- pcp
%      2 -- tmp
%      3 -- RO  runoff
%      4 -- Z index
%      5 -- X  PDSI
%      6 -- W  soil moisture (mid month)
%      7-16 SI  stress index, which is P-cPE, where c = .1 ,.2, ... 1.0
%
%   **Lomit is logical cv telling which tree sites to omit from use in the 
%     the modeling.  Lomit is a logical column vector of length nsites.  
%     An entry of 1 means omit that site from modeling.  A site might be omitted
%     because of unreliable data (e.g., chronology based on 1 core), or because
%     a previous run of reglag1 showed that the site has negligible climate signal.
%   **npred  is integer cv specifying how many predictors to allow to enter
%     npred is operative only when krule==3;  
% path\filename of .mat file with text information on tree-ring sites 
%     This file is generated by sitstr01.m, as described in \...steps.doc
% path\filename of .txt file with text info on sites (needed for elevations)
% path\filename of .dat file with lat/long of sites
%   
%----  OTHER NOTES 
%
% reglag1.m written specifically for the Sacramento River reconstruction. 
% reglag1.m  % has several goals.  When kmode==1, evaluates which of several types of water
% balance variables are best related to tree-ring series. When kmode==2, filters and
% scales the tree-ring series proportional to the strength of their moisture signal
% as a preliminary step to running PC regression using recflow1.m
%
% Deals with any of 16 different water-balance variables. Sixteen because
% that is how many 'output' water-balance variables my Palmer Drought Index program 
% generates.
%
% If kmode==2, models are validated
% in two ways:  split-sample and leave-n-out cross-validation.
%** Split Sample.  First find the subset of years in yrgo to yrsp with valid overlap
% of climate data and lagged tree-ring data.  Then split that subset into two halves.
% If total number of years is odd, second half if largest.
% If total number of years is even, first half and second half are equal in length.
% Calibrateson second half and verifies on first. Then calibrate on first half and
% verifies on second.  Fit entire data period for final "full period" model.
%** For cross-validation, uses leave-n-out strategy, where n is an odd number 
% depending on the maximum positive or negative lags used in the regression model.
% Number to 'leave out' is computed as 1 + 4*nlag, where nlag is the maximum allowable
% positive or negative lags.  If no lags are used, defaults to leave-1-out.  If 
% maximum lag is 2, then it's leave-9-out.
%
% Missing data handling.  The period yrgo-yrsp is a stipulated 
% longest possible modeling period. No years outside [yrgo yrsp] are considered
% for modeling.  Each tree-ring series might cover different years.  Each site-centered
% climate series might also cover different years.  The years within
% [yrgo yrsp] with both types of data define the full period for modeling.  Allowance is
% made for the need for 'endpoint' data for lagged predictors.
%
% Assumed form of tree-ring data.  I assume that arwhite2.m has been used to 
% whiten the tree-ring indices.  The tree-ring .mat file has the whitened indices as
% well as information on the AR whitening models -- model order and pct variance
% due modeled autocorrelation.
%
% Assumed naming of water-balance files.  The water-balance data is assumed to
% have been previously stored in .mat files, one file for tree-ring site. The 
% names are assumed to be wbout1.mat, wbout2.mat, ... because that is the convention
% used by my Palmer Index program.
%
% Use of npred.  Here is a typical sequence.  Say there are 12 sites. Set npred as
% ones(12,1), and run the function, saving results.  Then change npred to a vector
% of twos and repeat. Same for vector of 3s, 4s, and 5s, which would cover all 
% potential predictors in, say, a +-lag2 model. The examine outputs to determine
% which model is best for each site.  One way I have done this is to require:
%
%  R-squared >0.1
%  RE>0
%  REi+1 - Rei >0 (RE must increase with each variable entered)
%  (RSq - RE)>.05 (not too much drop off from calib to valid stats)
%------------- STEPS
%
% * READ INPUT AND INITIALIZE
% * SET DATA APPLICABLE TO ALL TREE-RING SITES
% * BUILD LAGGED PREDICTOR MATRIX
% * LOOP OVER TREE-RING SITES
%   * COMPUTE PERIOD OF VALID DATA
%   * IF KMODE==1, LOOP OVER 16 WATER-BALANCE VARIABLES
%   * IF KMODE==2, PULL THE DESIRED WATER-BALANCE VARIABLE
%		* FULL PERIOD MODEL
%		* OPTIONAL R-SQUARED COMPUTATION FOR FULL-PERIOD MODEL RE-RIT TO
%						COMMON PERIOD
%		* LONG TERM RECONSTRUCTION
%		* LEAVE-N-OUT CROSS-VALIDATION
%		* SPLIT-SAMPLE VALIDATION
% * CLEAN UP SOME VARIABLES
% * OPTIONALLY RE-SCALE FILTERED TREE-RING SERIES 
% * BUILD OUTPUT TABLES
%
%
%************************* TABLE SUMMARY ***********************
%
% General.
% Tables include entries for only those sites with Lomit~=1.
%
% TABLE1 -- tree-ring site names, .crn files and coordinates
%
%     Sequential Site no. in table (and optionally in original matrix)
%		.crn file name
%     Site name (first 8 chars), state, species
%		lat, long, elev (m)
%		start, end yrs of chron (before whitening & any other adjustment)
%
% TABLE 2 -- reconstruction summary
%
%     Site no (sequential in Table)
%		site name, state, species, and seq number (in original file)
%     First and last year of full reconstruction period
%     Lag structure of model
%		R-sqd, full-period regression model
%     RE, cross-validation
%
%     footnotes:  
%      reconstructed variable
%      Range of overall-F levels and their p-values
%
% TABLE 3 -- prewhitening models
%
%     Site no.
%     Site name, state species
%     ARmodel
%     AR var expld
%
% TABLE 4 -- regression structure and coefficients
%
%    Site #
%    Site name8, state, species code
%    Lag structure
%	  Full-period model reg coefs: const ... -2, -1 0 1 2 ...
%
%
% TABLE 5 -- RE summary
%
%    Site #
%    Site state, species code
%    R-squared, full-period reconstruction
%    RE  , cross-validation
%    R-Sq, calibration, for second-half calibration
%    RE  , validation, for second-half calibration
%    R-Sq, calibration, for first-half calibration
%    RE  , validation, for first-half calibration
%
%
% And this -- for RMSE summary TABLE 6
%
%    Site #
%    Site name
%    RMSE, full-period reconstruction model
%     "  , cross-validation
%     "  , calibration, for second-half calibration
%     "  , validation, for second-half calibration
%     "  , calibration, for first-half calibration
%     "  , validation, for first-half calibration
%
%
% TABLE 1------- for exploratory mode only
%
%    Site # (with original site number before Lomit screening in parentheses)
%    Sitename8, state, species
%    First and last years of full calibration period
%    R-Squared, PPT
%    R-Squared, Z index
%    R-Squared, PDSI
%    R-Squared, Temperature
%    R-Squared, average soil moisture
%  
%***********************************
% TEXT SUMMARY OF ANALYSIS (example)-- TXTSUM
% date/time of analysis
% Run of reglag1.m on ??.ctl
% Tree ring indices in ??
% Water Balance Variables in ??
% Exploratory mode
% Scaling called for
% Full-model period inclusive limits 1890 1995
% Common period for scaling:  NaN NaN
% Number of tree-ring sites in input matrix
% Number of tree-ring sites processed   = 
% Number of negative and positive lags on tree-rings 2,2 
%  
%

tic; % start timer to check how long to run this function

% A NaN variable to fill 3-dim arrays
a=NaN;

%---------- READ INPUT CONTROLS
[file1,path1]=uigetfile('*.ctl','Input control for reglag1.m');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

% Number of tree ring sites
nsites = str2num(fgetl(fid1)); % number of tree-ring sites

% File of the tree-ring data
pf2 = fgetl(fid1);
eval(['load ' pf2]);
if ~exist('X','var');
   error('X does not exist');
end
[mX,nX]=size(X);
if nsites ~= (nX-1);
   error('Specified nsites inconsistent with col size of X');
end
% store year column for X, and remove year col from X
yrX  = X(:,1);
yrsX=[min(yrX) max(yrX)];
X(:,1)=[];

% Scale tree-ring indices to usual "mean of 1" rather than 1000
X = X/1000;

% Path to the watbal output files
path3=fgetl(fid1);

% Season code specifying name of matrix of seasonalize water balance output
seascode=fgetl(fid1);

%----- Option for cutting off entry of predictors in regression
ctemp=fgetl(fid1); % krule (either 1,2 or 3);
krule=str2num(ctemp);
if krule~=1 & krule~=2 & krule~=3;
   error('krule must be 1,2 or 3');
end
ctemp=fgetl(fid1);
deltarsq = str2num(ctemp);
if deltarsq<=0 or deltarsq>=1.0;
   error('deltarsq must be greater than 0 and less than 1');
end

%-start end year of most inclusive period to be considered viable for modeling
ctemp = fgetl(fid1);
% get rid of leading and trailing blanks
ctemp=deblank(ctemp);
ctemp=deblank(fliplr(ctemp));
ctemp=fliplr(ctemp);
% take start as year ended by first blank
yrgo = str2num(strtok(ctemp));
% cut off chars thru first space
i1 = find(isspace(ctemp));
ctemp(1:i1(1))=[];
% get rid of leading blanks
ctemp=fliplr(deblank(fliplr(ctemp)));
yrsp = str2num(ctemp);
nyrgosp = yrsp-yrgo+1;  % number of years in uniform calibration period

%-start end year of period for computing means and standard devs
ctemp = fgetl(fid1);
% get rid of leading and trailing blanks
ctemp=deblank(ctemp);
ctemp=deblank(fliplr(ctemp));
ctemp=fliplr(ctemp);
% take start as year ended by first blank
yrstango = str2num(strtok(ctemp));
% cut off chars thru first space
i1 = find(isspace(ctemp));
ctemp(1:i1(1))=[];
% get rid of leading blanks
ctemp=fliplr(deblank(fliplr(ctemp)));
yrstansp = str2num(ctemp);

% Check whether scaling of filtered series by R-squared desired
if isnan(yrstango) | isnan(yrstansp);
   kscale = 0;
   ncomm=NaN;
else;
   kscale=1;
   if yrstango<yrgo | yrstansp>yrsp;
      error('yrstango and yrstansp must be within [yrgo yrsp]');
   end
   if (yrstansp-yrstango)<=20;
      error('Period defined by [yrstango yrstansp] must cover at least 21 yr');
   end
   ncomm = yrstansp-yrstango+1; % number of years in common period
end

% Number of positive and neg lags
nlag=str2num(fgetl(fid1));

% Mode (exploratory regression or final filtering)
kmode = str2num(fgetl(fid1));


% File with cv telling which clim variable to reconstruct, which sites to omit
% from regression analysis, and maximum number of predictors to allow in model
pf4 = fgetl(fid1);
eval(['load ' pf4]);
if ~exist('iclim','var');
   error([pf4 ' does not contain iclim']);
else
   if size(iclim,1) ~=1 | size(iclim,2)~=1;
      error(['iclim in ' pf4 ' must be a  scalar']);
   end
end
switch iclim;
case 1;
   climvar='Precip';
case 2;
   climvar='Temper';
case 3;
   climvar='RO';
case 4;
   climvar='Zindex';
case 5;
   climvar='PDSI';
case 6;
   climvar='Soilm'; % mid month soil moisture
case 7;
   climvar='P-0.1PE'; % stress index
case 8;
   climvar='P-0.2PE';
case 9;
   climvar='P-0.3PE';
case 10;
   climvar='P-0.4PE';
case 11;
   climvar='P-0.5PE';
case 12;
   climvar='P-0.6PE';
case 13;
   climvar='P-0.7PE';
case 14;
   climvar='P-0.8PE';
case 15;
   climvar='P-0.9PE';
case 16;
   climvar='P-1.0PE';
otherwise
end

% Lomit is column vector.  Ones mean omit that tree-ring site from analysis
if ~exist('Lomit','var');
   error([pf4 ' does not contain iclim']);
else
   if size(Lomit,1) ~=nsites | size(Lomit,2)~=1 | ~islogical(Lomit);
      error(['Lomit in ' pf4 ' must be a  logical cv of length nsites']);
   end
end

% npred is column vector telling maximum number of predictors to allow to enter
% for each site
if ~exist('npred','var');
   error([pf4 ' does not contain npred']);
else
   if size(npred,1) ~=nsites | size(Lomit,2)~=1 
      error(['npred in ' pf4 ' must be a  cv of length nsites']);
   end
end

 
% File with matrix of site string information
pf7 = fgetl(fid1);
eval(['load ' pf7]);
if ~exist('Snames','var');
   error([pf7 ' does not contain Snames']);
else
   if size(Snames,1) ~=nsites;
      error(['Snames in ' pf7 ' does not have ' int2str(nsites) ' rows'])
   end
end

% .txt file with site info, needed for elevation in meters
pf8 = fgetl(fid1);
eval(['elevm=load(pf8);']);
if size(Snames,1) ~= length(elevm);
   error(['number of elevations from ' pf8 ' not equal to number of sites in Snames ']);
end


% .dat file with long/lat of tree-ring sites
pf9 = fgetl(fid1);
eval(['lonlat=load(pf9);']);
if size(Snames,1) ~= size(lonlat,1);
   error(['no. of long/lats in ' pf9 ' not equal to number of sites in Snames ']);
end

% .dat file with start and end years of tree-ring chrons.  These unadjusted for
% 'ASS' years 
pf10 = fgetl(fid1);
eval(['yrgosp=load(pf10);']);
if size(Snames,1) ~= size(yrgosp,1);
   error(['no. of rows in ' pf10 ' not equal to number of sites in Snames ']);
end

%****************  SET DATA APPLICABLE TO ALL TREE-RING SITES

% -------------Allocate.  If final reconstruction mode, some variables will be
% column vectors, while if in exploratory mode they will be 16-column vectors.
% C4 will be either a 2-dim matrix or 3-dimensional array, again depending on
% kmode.  Note the "a" and "b" endings on variables.  "a" refers to model
% validated on first half and calibrated on second half.  "b" refers to model
% validated on second half and calibrated on first half. Neither "a" nor "b" means
% for the full-period model.  Likewise, C5 will hold
% 2SE values for the coef estimates.

if kmode==1;
   % Calibration R-squared
   C1 = repmat(NaN,nsites,16);  % full-period
      
   % Calibration overall F for regression eqn
   C2 = C1; % Calib F-level, overall
      
   % Calibration p-value for overall F
   C3 = C1; 
      
   % Regression coefficients
   C4 = repmat(NaN,[nsites 16 (2*nlag+2)]);
      
   % 2 standard devs of the coefs
   C5 = C4; 
      
   % Calibration period start and end years, and number of years in period
   YRSfull=repmat(NaN,nsites,2);  % first, last year for calibration
   Nfull =repmat(NaN,nsites,1); % number of years for calibration
   
   % Lag-structure matrix.  This is a coded string indicating which lags on the tree
   % rings are included in the final full-period regression model
   LAGSTR = repmat(blanks(7),nsites,1);

   
   % The following not used in exploratory analysis
   C1comm=[];
   V1a=[]; V2a=[]; V1b=[]; V2b=[]; V3a=[]; V3b=[];
   C1median=[];
   MAE=[]; RMSE=[]; RE=[];
   NP=[];
      
   
elseif kmode==2;
   
   % Calibration R-squared
   C1 = repmat(NaN,nsites,1);  % full-period
   C1a = C1; % validate on first half, calibrate on second
   C1b = C1; % validate on second half, calibrate on first
   
   % Calibration overall F for regression eqn
   C2 = C1; % Calib F-level, overall
   C2a = C2;  
   C2b = C2; 
   
   % Calibration p-value for overall F
   C3 = C1; 
   C3a = C3;
   C3b = C3; 
   
   % Regression coefficients
   C4 = repmat(NaN,nsites,(2*nlag+2));
   C4a = C4;
   C4b = C4; 
   
   % 2 standard devs of the coefs
   C5 = C4;
   C5a = C5; 
   C5b = C5; 
   
   % Rsquared for full-period model re-fit to common period.  The model's predictors
   % are set as those selected in calibrating the full-period model. Then the subset
   % of data for the common period is used to estimate the same coefficients again, and
   % to get the other regression statistics.
   C1comm = repmat(NaN,nsites,1);    
   
      
   %--------- Split sample statistics
   
   % root-mean-square error (rmse) of validation
   V1a =C1a; 
   V1b = C1a; 
   
   % reduction of error statistic (RE)
   V2a = C1a; 
   V2b = C1a; 
   
   % mean absolute error
   V3a = C1a; 
   V3b = C1a; 
   
   %-------- cross-validation statistics
   
   % median R-squared of cross-validation models for each site
   C1median = repmat(NaN,nsites,1);  

   %------  Other full-period regression statistics for reconstruction model
   
   YRSfull=repmat(NaN,nsites,2);  % first, last year for calibration
   Nfull =repmat(NaN,nsites,1); % number of years for calibration
   NP = Nfull;  % number of predictors (not counting constant) in final model
   
   
   %------- Some hybrid matrices for comparing calibration and validation stats
   MAE = repmat(NaN,nsites,2); % mean-absolute-errorfor full-period regression (col 1), 
   %     for cross-validation predictons (col 2)
   RMSE = repmat(NaN,nsites,2); % root-mean-square error for full-period regression (col 1), 
   %     for cross-validation predictons (col 2)
   RE = repmat(NaN,nsites,2); % reduction of error for full-period regression (col 1), 
   %     for cross-validation predictons (col 2).  For full-period regression, RE is
   %     identical to R-squared for regression
   
   % Lag-structure matrix.  This is a coded string indicating which lags on the tree
   % rings are included in the final full-period regression model
   LAGSTR = repmat(blanks(7),nsites,1);
   
   
end;  % of for kmode==

   

%************** BUILD LAGGED PREDICTOR MATRIX

[XL,yrXL] = lagyr3(X,yrsX,[0 nlag nlag]);
yrXL = (yrXL(1,1):yrXL(1,2))';
% note: yrXL not used in this pgm

% optional pointer to common period years in XL
if kmode==2 & kscale==1;
   LXLcomm = yrXL>=yrstango & yrXL<=yrstansp;
else
   LXLcomm=[];
end

% Store uniform calib pd years of lagged predictor matrix.  Will need this in both
% exploratory and reconstruction mode. 
Ltemp = yrX>=yrgo  & yrX<=yrsp;
U = XL(Ltemp,:);
yrU = (yrgo:yrsp)';

% If in exploratory regression mode, can get by with only the short data in XL and yr
if kmode==1;
   XL = U;
   yrXL = yrU;
end

% Size the matrix of lagged tree-ring data.  This size might have changed if you
% selected exploratory mode.
[m2,n2]=size(XL);

% A constant column vector of ones for use in regression modeling
onecv = ones(m2,1); 

% Storage for the filtered output tree-ring series, needed depending on kmode
if kmode==1;
   YL=[];
else;
   YL = repmat(NaN,m2,nsites);  
end

% Clear unneded variables
clear X Ltemp


%---- Prompt user for names of output files

if kmode==1;
   [file6,path6]=uiputfile('rlo*.mat','Store output from exploratory-mode here');
   pf6=[path6 file6];
   %fid6=fopen(pf6,'w');
   %fprintf(fid6,'%s',head6);
else;
   [file5,path5]=uiputfile('rlo*.mat','Store output from final-mode here');
   pf5=[path5 file5];
   
end

%---- Prompt user for output Table options
kseq = questdlg('Include original site sequence number in output tables?');
if strcmp(kseq,'Cancel');
   kseq='No';
end
% kseq=='Yes' means include it

%************** LOOP OVER TREE-RING SITES

clc
disp('STARTING LOOP OVER TREE-RING SITES');

% Have total of nsites, but previous info may have led you to mark some
% sites for omission from analysis (by Lomit==1).  Compute integer column 
% vector marking sites to analyze.  Also compute total number of sites you
% want to analyze. Do only sites for which Lomit~=1;
nslean = (1:nsites);  % vector for all sites
if any(Lomit);
   nslean = nslean(~Lomit);
end
nsites2 = length(nslean); % number of sites to analyze (with Lomit~=1)

% Loop over the to-be-used tree-ring sites
for ns=nslean; 
   disp(['Site # ' int2str(ns)]);
   
   % Recall that XL is matrix of lagged tree-ring indices for all sites. Make
   % a row vector indexing the columns of XL holding data for this tree-ring site
   mintemp = ns;
   maxtemp = ns + (nlag*2)*nsites;
   icol1 = mintemp:nsites:maxtemp;
   
   % Cleanup
   clear mintemp maxtemp;
   
	%****** COMPUTE PERIOD OF VALID DATA 

   % Compute first, last year for which the submatrix of current and lagged
   % tree rings for this site has no missing data
   Xtemp = XL(:,icol1);
   if length(icol1)>1;
      Ltemp = (~any(isnan(Xtemp')))';
   else;
      Ltemp = ~isnan(Xtemp);
   end
   yrgoodx = yrXL(Ltemp); % col vector of years
   yrsgoodx=[min(yrgoodx) max(yrgoodx)]; % first, last year
   
   % Get the water-balance file for this site
   file3=['wbout' int2str(ns)];
   pf3=[path3 file3];
   eval(['load ' pf3]);
   clear datout; % do not need this huge cell variable
   % Put the matrix of climate variables in W
   eval(['W = ' seascode ';']);  % pull the desired matrix of seasonalized time series
   
   % Get years for water balance matrix
   yrW = W(:,1);
   yrsW = [min(yrW) max(yrW)];
   W(:,1)=[]; % drop year col
   
   % Note that columns of W are now:
   %  vbl 1 -- pcp
   %   vbl 2 -- tmp
   %   vbl 3 -- RO  runoff
   %   vbl 4 -- Z index
   %   vbl 5 -- X  PDSI
   %   vbl 6 -- W  soil moisture (mid month)
   %   vbls 7-16 SI  stress index, which is P-cPE, where c = .1 ,.2, ... 1.0
   
   
   %****** IF KMODE==1, LOOP OVER 16 WATER-BALANCE VARIABLES
   
   if kmode ==1;
      for nv = 1:16; % Loop over water balance variables
         w = W(:,nv); % a water-balance time series
         
         % Get years for which w has valid water-balance data
         Ltemp =~isnan(w);
         yrgoodw = yrW(Ltemp);
         yrsgoodw = [min(yrgoodw) max(yrgoodw)];
         
         % Mark end yrs of period within [yrgo yrsp] for which have
         % valid water-balance data and and lagged tree-ring data
         % This is defined as the "full period" for modeling
         yrs3= [max([yrsgoodw(1) yrsgoodx(1) yrgo])  min([yrsgoodw(2) yrsgoodx(2) yrsp])];
         yr3 = (yrs3(1):yrs3(2))';
         nyrs3 = length(yr3);
          
                  
         % Make pointers to the three analysis periods in w and XL
         Lw3 = yrW>=yrs3(1) & yrW<=yrs3(2);
         LXL3 = yrXL>=yrs3(1) & yrXL<=yrs3(2);
                
          % Store calibration period years for this site
          YRSfull(ns,:)= yrs3;
          Nfull (ns) = nyrs3;
         
         %-------  Full-Period Analysis (the only analysis done for kmode==1)
         
         % Pull off the full-period data
         Xkey = XL(LXL3,:);
         wkey = w(Lw3);
         
         % Do the regression 
         if krule==1;
            [I2,I4,stats,c,e,what]=stepr2(Xkey,wkey,icol1);
         elseif krule==2;
            [I2,I4,stats,c,e,what]=stepr2(Xkey,wkey,icol1,deltarsq);
         else
            error('krule must be 1 or 2 for a mode-1 analysis');
         end
         
         % Lag structure
         s=blanks(7);
         if any(I4==1); s(1)='0'; end;
         if any(I4==2); s(2)='N'; s(3)='1'; end;
         if any(I4==3); s(2)='N'; s(4)='2'; end;
         if any(I4==4); s(5)='P'; s(6)='1'; end;
         if any(I4==5); s(5)='P'; s(7)='2'; end;
         iblank = findstr(s,' ');
         if ~isempty(iblank);
            s(iblank)=[];
            nadd = 7-length(s);
            s = [s blanks(nadd)];
         end
         LAGSTR(ns,:)=s;

         
         
         % Store results
         C1(ns,nv)=stats(1); % R-sqd for full-period model
         C2(ns,nv)=stats(3); % overall F for equation
         
         % P value for overall F   
         Fcompute=stats(3);
         Fp=1-fcdf(Fcompute,(length(c)-1),(nyrs3-length(c)));
         C3(ns,nv)=Fp;
      
         % Coefficients
         C4(ns,nv,1)=c(1);
         C4(ns,nv,I4+1)=c(2:length(c));
         
      end; % of loop nv over water balance variables
      
      
      %********* IF KMODE==2, PULL THE DESIRED WATER-BALANCE VARIABLE
      
      
   elseif kmode==2;  % Final filtering mode
      w = W(:,iclim); % the pre-selected water-balance time series
      
      % Get years for which w has valid data
      Ltemp =~isnan(w);
      yrgoodw = yrW(Ltemp);
      yrsgoodw = [min(yrgoodw) max(yrgoodw)];
      
      % Mark end yrs of period within [yrgo yrsp] with clim data and valid lagged
      % tree-ring data
      % This is the full period for modeling
      yrs3= [max([yrsgoodw(1) yrsgoodx(1) yrgo])  min([yrsgoodw(2)  yrsgoodx(2) yrsp])];
      yr3 = (yrs3(1):yrs3(2))';
      nyrs3 = length(yr3);
      
      % If scaling desired, modeling period must include [yrstango yrstansp]
      if kscale==1;
         if yrstango<yrs3(1) | yrstansp>yrs3(2);
            error(['   Model years yrs3 do not contain common period']);
         end
      end
      
      % Make year variables for first and second calib/valid periods
      nleft = rem(nyrs3,2);
      if nleft==0; % even number of years
         nyrs1=nyrs3/2;
         nyrs2=nyrs1;
         nyrs1 = ceil(nyrs3/2);
         nyrs2 = nyrs3-nyrs1;
         yr2 = yr3(1:nyrs2);
         yr1 = yr3((nyrs2+1):nyrs3);
         yrs2 = [min(yr2) max(yr2)];
         yrs1 = [min(yr1) max(yr1)];
      else; % odd number of years
         nyrs1 = ceil(nyrs3/2);
         nyrs2 = nyrs3-nyrs1;
         yr2 = yr3(1:nyrs2);
         yr1 = yr3((nyrs2+1):nyrs3);
         yrs2 = [min(yr2) max(yr2)];
         yrs1 = [min(yr1) max(yr1)];
      end
      
      % Make pointers to the three analysis periods in w and XL
      Lw3 = yrW>=yrs3(1) & yrW<=yrs3(2);
      LXL3 = yrXL>=yrs3(1) & yrXL<=yrs3(2);
      Lw2 = yrW>=yrs2(1) & yrW<=yrs2(2);
      LXL2 = yrXL>=yrs2(1) & yrXL<=yrs2(2);
      Lw1 = yrW>=yrs1(1) & yrW<=yrs1(2);
      LXL1 = yrXL>=yrs1(1) & yrXL<=yrs1(2);
      % common period
      Lwcomm = yrW>=yrstango & yrW<=yrstansp;
            
      % Pull off the full-period data
      Xkey = XL(LXL3,:);
      wkey = w(Lw3);
      
      % Store calibration period years for this site
      YRSfull(ns,:)= yrs3;
      Nfull (ns) = nyrs3;
      
      
      %********* FULL PERIOD MODEL
      
      if ns==63; 
         disp('site 63');
      end
      
      
      if krule==1;
         [I2,I4,stats,c,e,what]=stepr2(Xkey,wkey,icol1); % do regression
      elseif krule==2;
         [I2,I4,stats,c,e,what]=stepr2(Xkey,wkey,icol1,deltarsq);
      elseif krule==3;
         [I2,I4,stats,c,e,what]=stepr2(Xkey,wkey,icol1,deltarsq,npred(ns));
      else
         error('krule must be 1,2 or 3');
      end
                  
      % Store results
      C1(ns)=stats(1); % R-sqd for full-period model
      C2(ns)=stats(3); % overall F for equation
      
      % P value for overall F   
      Fcompute=stats(3);
	   Fp=1-fcdf(Fcompute,(length(c)-1),(nyrs3-length(c)));
      C3(ns)=Fp;
      
      % Coefficients
      C4(ns,1)=c(1);
      C4(ns,I4+1)=(c(2:length(c)))';
      
      % Lag structure
      s=blanks(7);
      if any(I4==1); s(1)='0'; end;
      if any(I4==2); s(2)='N'; s(3)='1'; end;
      if any(I4==3); s(2)='N'; s(4)='2'; end;
      if any(I4==4); s(5)='P'; s(6)='1'; end;
      if any(I4==5); s(5)='P'; s(7)='2'; end;
      iblank = findstr(s,' ');
      if ~isempty(iblank);
         s(iblank)=[];
         nadd = 7-length(s);
         s = [s blanks(nadd)];
      end
      LAGSTR(ns,:)=s;
               
      %------ A few other full-model calibration statistics
      xfull = [ones(nyrs3,1) Xkey(:,I2)]; % predictor matrix
      yhfull = xfull * c; % predictions
      [dopey,MAE(ns,1),RMSE(ns,1),RE(ns,1)]=cvstat3(wkey,yhfull,repmat(mean(wkey),nyrs3,1));
      clear xfull yhfull dopey;
      NP(ns)=length(I2); % store number of predictors in model
      
      %******** OPTIONAL R-SQUARED COMPUTATION FOR FULL-PERIOD MODEL RE-FIT TO COMMON PERIOD
      if kscale==1;
         % recall that LXLcomm is pointer to common period in lagged tree-ring matrix
         % and Lwcomm  to common period in climate variable
         
         % Get common-period data
         xcomm = [ones(ncomm,1) XL(LXLcomm,I2)]; % predictor
         wcomm = w(Lwcomm);
         
         % Regression
         [Bdum,BINTdum,Rdum,RINTdum,STATSdum] = REGRESS(wcomm,xcomm,.01);
         C1comm(ns)=STATSdum(1);% store R-squared
         clear Bdum BINTdum Rdum RINTdum STATSdum;
      end
            
      
         
      %************ LONG TERM RECONSTRUCTION
      XTEMP = [onecv XL(:,I2)]; % predictor matrix
      y =  XTEMP*c;  % reconstruction == filtered series
      YL(:,ns)=y; % store
      
      disp(['  Full-model and regression done']);
      
      
      %******  LEAVE-N-OUT CROSS-VALIDATION
      
      whcv = repmat(NaN,nyrs3,1); % to store cross-validation estimates
      wcvmean = repmat(NaN,nyrs3,1); % to store cross-validation calibration
      % subset means of predictand
      crsq = repmat(NaN,nyrs3,1);  % to store R-squared of calib for each
      % cross-validatio model
      IX=crospull(nyrs3,nlag,nlag); % logical pointer to subsets of years
      
      % Loop over the cross-validation models for this site
      for ncv = 1:nyrs3;
         % Get data for the subperiod years
         Lcv = IX(:,ncv);
         nyrcv = sum(Lcv); % number of valid years for this calibration
         xcv = Xkey(Lcv,:);
         wcv = wkey(Lcv); % predictand
         wcvmean(ncv) = mean(wcv);  % calib period mean of obs data for sub period
         
         % Estimate regression coefs and get regr stats
         [c,stats]=sos(wcv,xcv(:,I2));
         
         % Store R-squared of calibration for each cross-validation model
         crsq(ncv) = stats(1);
         
         % Apply regression to get predictions, for full calibration years
         xtemp = [ones(nyrs3,1)  Xkey(:,I2)]; % predictor matrix
         ypredcv = xtemp * c; % predictions
         whcv(ncv) = ypredcv(ncv); % cross-validation estimate
      end
               
      % Compute cross-validation statistics
      [r2,mae,rmse,re]=cvstat3(wkey,whcv,wcvmean);
      MAE(ns,2)=mae;
      RMSE(ns,2)=rmse;
      RE(ns,2)=re;
      C1median(ns) = median(crsq);
      
      % Cleanup
      clear r2 mae rmse re
      
      disp(['  Cross-validation done']);
      
      
      
      %********** SPLIT-SAMPLE VALIDATION
      
      
      %---- Calibrate on second half, validate on first half
      
      % Pull calibration period data
      Xkey = XL(LXL1,:);
      wkey = w(Lw1);
      wmeancal = mean(wkey); % calibration-period mean of predictand
      
      % Estimate regression coefs and get regr stats
      [c,stats]=sos(wkey,Xkey(:,I2));
      
      % Store calibration statistics
      C1a(ns)=stats(1); % R-sqd 
      C2a(ns)=stats(3); % overall F for equation
      
      % P value for overall F   
      Fcompute=stats(3);
	   Fp=1-fcdf(Fcompute,(length(c)-1),(nyrs1-length(c)));
      C3a(ns)=Fp;
      
      % Regression coeffs
      C4a(ns,1)=c(1);
      C4a(ns,I4+1)=(c(2:length(c)))';
      
      % Pull validation data
      Xkey = XL(LXL2,:);
      wkey = w(Lw2);
      
      % Apply regression to get predictions for validation period
      xtemp = [ones(nyrs2,1)  Xkey(:,I2)]; % predictor matrix
      whvalid = xtemp * c; % predictions

      % Compute validation statistics
      [r2,V3a(ns),V1a(ns),V2a(ns)]=cvstat3(wkey,whvalid,repmat(wmeancal,nyrs2,1));
      clear r2
                        
      
      
      %---- Calibrate on first half, validate on second half
      
      % Pull calibration period data
      Xkey = XL(LXL2,:);
      wkey = w(Lw2);
      wmeancal = mean(wkey); % calibration-period mean of predictand
      
      % Estimate regression coefs and get regr stats
      [c,stats]=sos(wkey,Xkey(:,I2));
      
      % Store calibration statistics
      C1b(ns)=stats(1); % R-sqd 
      C2b(ns)=stats(3); % overall F for equation
      
      % P value for overall F   
      Fcompute=stats(3);
	   Fp=1-fcdf(Fcompute,(length(c)-1),(nyrs2-length(c)));
      C3b(ns)=Fp;
      
      % Regression coeffs
      C4b(ns,1)=c(1);
      C4b(ns,I4+1)=(c(2:length(c)))';
      
      % Pull validation data
      Xkey = XL(LXL1,:);
      wkey = w(Lw1);
      
      % Apply regression to get predictions for validation period
      xtemp = [ones(nyrs1,1)  Xkey(:,I2)]; % predictor matrix
      whvalid = xtemp * c; % predictions

      % Compute validation statistics
      [r2,V3b(ns),V1b(ns),V2b(ns)]=cvstat3(wkey,whvalid,repmat(wmeancal,nyrs1,1));
      clear r2
      
      disp(['  Split-sample validation done']);

   end ; %  of else kmode==2
   
   
end; % of loop ns over tree ring sites


if kmode==2;
%********* OPTIONALLY RE-SCALE FILTERED TREE-RING SERIES 
   if kscale==1;
      Z = YL(LXLcomm,:);
      [mz,nz]=size(Z);
      zmean = nanmean(Z);
      zsd = nanstd(Z);
      
      % Convert YL to zscores
      Y2 = (YL -(repmat(zmean,m2,1))) ./  (repmat(zsd,m2,1));
      
      % Scale full reconstruction by R-squared of regression for common period
      fact1=C1comm';
      fact1=sqrt(fact1); % Ensures that variance of output series equals regression Rsquared
      fact1 = repmat(fact1,m2,1);
      Y2 = Y2 .*  fact1;
   else;  % scaling not desired
      Y2=YL;
   end
   
   % Lop off rows corresponding to sites not analyzed (Lomit==1)
   LAGSTR=LAGSTR(~Lomit,:);
   
   
   % Put year column on storage matrix
   X2=[yrX Y2];
elseif kmode==1;
   % No action needed because no filtered or filtered, scaled series needed   
end

%------- TEXT SUMMARY OUTPUT

% Time to run function (up to this point)
timeoff=toc;
strtime = num2str(timeoff,'%8.1f seconds');

% Date
strdate =['Run of reglag1.m on ' datestr(now)];

% Entry of variables
if krule==1;
   strrule = 'krule=1:  adjusted R-squared must increase';
elseif krule==2;
   strrule =['krule=2: requires R-squared to increase by at least ' num2str(deltarsq)];
else;
   strrule =['krule=3: adj Rsq must increase and max num of predictors constrained by npred'];
end
   

str1='SUMMARY OF ANALYSIS';
str1=char(str1,' ');
str1=char(str1,strdate);
str1=char(str1,['   Time to run analysis = ' strtime]);
str1=char(str1,['Control info from ' pf1]);
str1=char(str1,['Tree ring data from ' pf2]);
str1=char(str1,['Path to input water balance data = ' path3]);
if kmode==2;
   str1=char(str1,['   Selected predictand = ' climvar]);
end
str1=char(str1,['   Seasonal grouping of predictand = ' seascode]);
str1=char(str1,['Lomit, iclim, npred from ' pf4]);
if kmode==1;
   str1=char(str1,['Exploratory output to ' pf6]);
else;
   str1=char(str1,['Final output to ' pf5]);
end
str1=char(str1,['Site names from ' pf7]);
str1=char(str1,['Elevations from ' pf8]);
str1=char(str1,['longitudes, lats of sites from ' pf9]);
str1=char(str1,['Start, end year of chrons from ' pf10]);
str1=char(str1,' ');
str1=char(str1,['Mode for analysis, kmode = ' num2str(kmode)]);
str1=char(str1,['Number of tree-ring sites in ' pf2 ' = ' int2str(nsites)]);
str1=char(str1,['Number analyzed (Lomit~=1) = ' int2str(nsites2)]);
str1=char(str1,['Max allowable neg and pos lags in regression = ' int2str(nlag)]);
str1=char(str1,['Rule for entry of another predictor in regression']);
str1=char(str1,['   ' strrule]);
if kscale==1;
   str1=char(str1,'Filtered output tree rings scaled by regression Rsqd on period');
   str1=char(str1,['   ' num2str([yrstango yrstansp])]);
else;
   if kmode==2;
      str1=char(str1,'No scaling of filtered tree-ring output desired');
   else; % exploratory mode
      str1=char(str1,'No filtered or scaled, filtered series produced in this mode');
   end;
end;

txtsum=str1;



%**********  RESTRICTED OUTPUT FOR EXPLORATORY MODE ONLY
if kmode==1;
   
   
   % TABLE 1------- for exploratory mode only
      
   head1a='No.    Site        Spec Period     PPT   Z     PDSI  Temp  Soilm';    
   head1b='No.    Site       Spec Period     PPT   Z     PDSI  Temp  Soilm';   
   if strcmp(kseq,'Yes'); % want original site sequence numbers in table also
      foot1 = 'site number -- sequential in table, (and original in input matrix)';
   else;
      foot1='sequential number of site in this table';
   end
   foot1 = char(foot1,'site name (first 8 chars) and state');
   foot1 = char(foot1,'species code');
   foot1 = char(foot1,'calibration period for regression');
   foot1 = char(foot1,'R-squared for regression with various predictand variables: ');
   foot1=char(foot1,'   PPT -- precipitation');
   foot1=char(foot1,'   Z -- Palmer Z-index');
   foot1=char(foot1,'   PDSI -- Palmer drought severity index');
   foot1=char(foot1,'   Temp -- average temperature');
   foot1=char(foot1,'   Soilm -- soil mositure at mid-month');
   
   bk1=repmat(' ',nsites2,1);
   bk2=repmat('  ',nsites2,1);
   bk5=repmat('     ',nsites2,1);

   
   % R-squared for regression models
   Cgroup = C1(:,[1 4 5 2 6]);
   strRSQ = num2str(Cgroup(~Lomit,:),'%6.2f'); % R-squared for 5 water balance variables
   
   strnum1 = num2str((1:nsites2)'); % site number, sequential in table
   strname8 = Snames(~Lomit,17:24);  % first 8 chars  of site name
   strstsp = Snames(~Lomit,40:46); % state and species code
   Lparen = repmat('(',nsites2,1);
   strnum2 = Snames(~Lomit,1:2); % site number in orignal tree-ring file
   Rparen = repmat(')',nsites2,1);
   str5 = [num2str(YRSfull(~Lomit,1),'%4.0f') repmat('-',nsites2,1) num2str(YRSfull(~Lomit,2),'%4.0f')];
   
   strall1a = [strnum1 Lparen strnum2 Rparen bk2 strname8 bk1 strstsp bk1];
   strall1b = [strnum1 bk5 strname8 bk1 strstsp bk1];
   strall1c = [str5 bk2 strRSQ];
   
   if strcmp(kseq,'Yes'); % want original site sequence numbers in table also
      strall=[strall1a strall1c];
      head1=head1a;
   else;
      strall = [strall1b strall1c];
      head1=head1b;
   end
   table1 = char(head1,' ',strall,' ',foot1);
   
   %-- List of saved variables 
   txtvars = 'Lomit -- logical pointer to tree-ring sites to omit from analysis';
   txtvars=char(txtvars,'nslean -- index to tree-ring sites analyzed. Index refers');
   txtvars=char(txtvars,'  to rows in xy coordinate array and other arrays, and can');
   txtvars=char(txtvars,'  be used to extract text information, long, latitude, and ');
   txtvars=char(txtvars,'  other information for the subset of sites analyzed');
   txtvars =char(txtvars,'C1 -- R-squared values for regression. Cols are for predictands: ');
   txtvars=char(txtvars,'  1 PPT');
   txtvars=char(txtvars,'  2 Temp');
   txtvars=char(txtvars,'  3 Runoff (model)');
   txtvars=char(txtvars,'  4 Palmer Z index');
   txtvars=char(txtvars,'  5 Palmer drought severity index');
   txtvars=char(txtvars,'  6 soil moisture at mid month');
   txtvars=char(txtvars,'  7-16 stress index defined as P-x*PE, where x - .1,.2,... 1.0');
   txtvars=char(txtvars,'C2 -- Overall F-level for regressions, in order as in C1');
   txtvars=char(txtvars,'C3 -- p-values for overall F-levels');
   txtvars=char(txtvars,'C4 regression coefficients for models, as follows:');
   txtvars=char(txtvars,'   constant, lag 0, lag t-1, lag t-2... lag t+1, lag t+2,....');
   txtvars=char(txtvars,'   C4 is 3-d, with dimensions for sites, variables and coefs');
   txtvars=char(txtvars,'   You can extract 2-d matrices of coefficients on all sites for');
   txtvars=char(txtvars,'   a specific variable or on all 16 variables for a specific site');
   txtvars=char(txtvars,'   by sqeezing the 3-d array as in following example: ');
   txtvars=char(txtvars,'      PDSI (variable 4) for all sites:');
   txtvars=char(txtvars,'         D4=C4(:,4,:);  E4=squeeze(D4);');
   txtvars=char(txtvars,'      Site 5,for all 16 variables');
   txtvars=char(txtvars,'         D4=C4(5,:,:);  E4=squeeze(D4);');
   txtvars=char(txtvars,'YRSfull first and last year of calibration periods');
   txtvars=char(txtvars,'txtsum text summary of settings for analysis');
   txtvars=char(txtvars,'txtvars  definitions of saved variables');
   txtvars=char(txtvars,'table1 ... -- summary table ');
   
   C1=C1(~Lomit,:);
   C2=C2(~Lomit,:);
   C3=C3(~Lomit,:);
   C4=C4(~Lomit,:,:);
   YRSfull=YRSfull(~Lomit,:);
      
   saveset1 = 'nslean Lomit YRSfull C1 C2 C3 C4 table1 txtvars txtsum';
   eval(['save ' pf6  ' ' saveset1 ' ;']);
   return;
end
   
%********** BUILD OUTPUT TABLES

% --- TABLE 2
%
%     Number (sequential)
%		site name, state, species, and seq number of site in original set
%     First and last year of full reconstruction period
%     Lag structure of model
%		R-sqd, full-period regression model
%     RE, cross-validation
%
%     footnotes:  
%      reconstructed variable
%      Range of overall-F levels and their p-values
head1a=' N  NAME         SPEC   PERIOD        LAGS    RSQ   RE';
head1b='N   NAME        SPEC   PERIOD        LAGS    RSQ   RE';

bk1 = repmat(' ',nsites2,1);
bk3= [bk1 bk1 bk1];

strnum1 = num2str((1:nsites2)'); % site number, sequential in table
strname = Snames(~Lomit,17:38);  % site name
strname8 = Snames(~Lomit,17:24);  % first 8 chars  of site name
strgo1 = num2str(yrgosp(~Lomit,1));

strstsp = Snames(~Lomit,40:46); % state and species code
Lparen = repmat('(',nsites2,1);
strnum2 = Snames(~Lomit,1:2); % site number in orignal tree-ring file
Rparen = repmat(')',nsites2,1);
str5 = [num2str(YRSfull(~Lomit,1),'%4.0f') repmat('-',nsites2,1) num2str(YRSfull(~Lomit,2),'%4.0f')];
str6 = LAGSTR;
str7 = num2str(RE(~Lomit,:),'%6.2f');
str8 = Snames(~Lomit,4:15);  % .crn file name

strall1a = [strnum1  Lparen strnum2 Rparen strname8 bk1 strstsp];
strall1b = [strnum1 bk1 bk1 bk1 strname8 bk1 strstsp];
strall2 = [bk1 str5 Lparen strgo1 Rparen bk1 str6 bk1 str7]; 

if strcmp(kseq,'Yes'); % want original site sequence numbers in table also
   strall=[strall1a strall2];
   head1=head1a;
else;
   strall = [strall1b strall2];
   head1=head1b;
end

foot1 = '1No.';
foot1 = char(foot1,'2site name, with state and species code');
foot1 = char(foot1,['xxsequence number among ' int2str(nsites) ' sites']);
foot1 = char(foot1,'3calibration period, with start year of chronology in parentheses');
foot1 = char(foot1,'4lags on predictors (tree-ring indices), coded as follows:');
foot1 = char(foot1,'  0= lag-zero, tree-ring index in same year as predictand');
foot1 = char(foot1,'  N1 = negative lag 1, index in year t-1  as predictor of climate year t');
foot1 = char(foot1,'  N2 = negative lag 2, index in year t-2 ...');
foot1 = char(foot1,'  N12 = negative lags 1 and 2, ...');
foot1 = char(foot1,'  P1, P2, P12 = as for N1,N2, N12, except tree-ring index lagged');
foot1 = char(foot1,'      forward from predictand');
foot1 = char(foot1,'5RSQ = R-squared for regression model');
foot1 = char(foot1,['6RE = reduction-of-error statistic from leave-' int2str(4*nlag+1) '-out']);
foot1 = char(foot1,'      cross-validation of regression model');
table2 = char(head1,strall,' ',foot1);


% --- TABLE1

head2a=' N    NAME                      SPEC  FILE         LONG      LAT   EL(M)  PERIOD';
head2b='N   NAME                      SPEC  FILE         LONG     LAT    EL(M)  PERIOD';

foot1='Site No. (with site number in original matrix in parens)';
foot1=char(foot1,'Site name and state');
foot1=char(foot1,'Species code');
foot1=char(foot1,'Longitude in decimal degrees');
foot1=char(foot1,'Latitude in decimal degrees');
foot1=char(foot1,'Elevation of site (meters)');
foot1=char(foot1,'First and last year of tree-ring chronology (before AR filtering)');

strfile = Snames(~Lomit,4:15); % .crn filenames
strxy = num2str(lonlat(~Lomit,:),'%7.2f %5.2f'); % long, lat
strgosp = num2str(yrgosp(~Lomit,:),'%5.0f'); % start, ent year of tree-ring chrons
strelm = num2str(elevm(~Lomit,:),'%5.0f'); % elev in meters

strall1a = [strnum1 Lparen strnum2 Rparen bk1 strname bk1 strstsp bk1 strfile bk1];
strall1b = [strnum1 bk3 strname bk1 strstsp bk1 strfile bk1];
strall2 = [strxy bk1 strelm bk1 strgosp bk1];

if strcmp(kseq,'Yes'); % want original site sequence numbers in table also
   strall=[strall1a strall2];
   head2=head2a;
else;
   strall = [strall1b strall2];
   head2=head2b;
end

table1 = char(head2,strall,' ',foot1);


%-----------  TABLE 3 (AR modeling)

head3a='N     SITE        SPEC Q VRAT';
head3b='N   SITE        SPEC Q VRAT';


strarord = num2str(arorder(~Lomit)); % order of AR prewhitening model
persist=1-arvrat; % % variance due modeled persistence
strvrat = num2str(persist(~Lomit),'%6.2f');

strall3a = [strnum1 Lparen strnum2 Rparen bk1 strname8 bk1 strstsp bk1];
strall3b = [strnum1 bk3 strname8 bk1 strstsp bk1];
strall3c = [strarord bk1 strvrat];

if strcmp(kseq,'Yes'); % want original site sequence numbers in table also
   strall=[strall3a strall3c];
   head3=head3a;
else;
   strall = [strall3b strall3c];
   head3=head3b;
end

foot1 = 'Order of AR prewhitening model';
foot1 = char(foot1,'% variance due modeled autocorrelation'); 

table3=char(head3,strall ,' ',foot1);


%  TABLE 4 -- lag structure of models

head4a='N     SITE        SPEC LAGS     const     0       -1     -2       +1      +2';    
head4b='N   SITE        SPEC LAGS     const     0       -1     -2       +1      +2';
foot1 = 'regression constant and coefficients on indicated lags';


strcoef = num2str(C4(~Lomit,:),'%8.3f');

strall4a = [strnum1 Lparen strnum2 Rparen bk1 strname8 bk1 strstsp bk1];
strall4b = [strnum1 bk3 strname8 bk1 strstsp bk1];
strall4c = [LAGSTR  bk1 strcoef];

if strcmp(kseq,'Yes'); % want original site sequence numbers in table also
   strall=[strall4a strall4c];
   head4=head4a;
else;
   strall = [strall4b strall4c];
   head4=head4b;
end

table4 = char(head4,' ',strall,' ',foot1);



%------ TABLE 5 - for RE summary, with comparison of cross-validation and split sample  results
%
%    Site #
%    Sitename8,, state, species code
%    R-squared, full-period reconstruction
%    RE  , cross-validation
%    R-Sq, calibration, for second-half calibration
%    RE  , validation, for second-half calibration
%    R-Sq, calibration, for first-half calibration
%    RE  , validation, for first-half calibration

head5a='N     SITE        SPEC RSQf  REf   RSQa  REa   RSQb  REb';     
head5b=head5a;
foot1 = 'R-squared for full-period calibration, with RE from cross-validation';
foot1 = char(foot1,'R-squared for calibration on late, validatio on early');
foot1 = char(foot1,'R-squared for calibration on early, validatio on late');


REtemp= [RE C1a V2a C1b V2b];
REcomp = REtemp;
strRE = num2str(REtemp(~Lomit,:),'%6.2f');

strall5a = [strnum1 Lparen strnum2 Rparen bk1 strname8 bk1 strstsp bk1];
strall5b = [strnum1 bk3 strname8 bk1 strstsp bk1];
strall5c = [strRE];

if strcmp(kseq,'Yes'); % want original site sequence numbers in table also
   strall=[strall5a strall5c];
   head5=head5a;
else;
   strall = [strall5b strall5c];
   head5=head5b;
end
clear REtemp

table5 = char(head5,' ',strall,' ',foot1);


%----------  SAVE OUTPUT FILES

if kmode==1;
   txt = 'C1 -- R-squared';
   txt=char(txt,...
      'C2 -- overall F',...
      'C3 -- p-value for overall F');
   eval(['save ' pf6  ' C1 C2 C3 txt ;']);
else;
   txtvars = 'Lomit -- logical pointer to tree-ring sites to omit from analysis';
   txtvars=char(txtvars,'nslean -- index to tree-ring sites analyzed. Index refers');
   txtvars=char(txtvars,'  to rows in xy coordinate array and other arrays, and can');
   txtvars=char(txtvars,'  be used to extract text information, long, latitude, and ');
   txtvars=char(txtvars,'  other information for the subset of sites analyzed');
   
   txtvars=char(txtvars,'NP number of predictor variables selected for each site');
   txtvars = char(txtvars,'X2 tsm of filtered and possibly scaled tree-ring series');
   txtvars =char(txtvars,'REcomp -- comparison matrix of Rsquared and RE values: ');
   txtvars=char(txtvars,'  1 Rsqd for full-period final reconstruction models');
   txtvars=char(txtvars,'  2 RE for cross-validation of full-period model');
   txtvars=char(txtvars,'  3 Rsqd for calibrate on late half');
   txtvars=char(txtvars,'  4 RE for validate on early half');
   txtvars=char(txtvars,'  5 Rsqd for calibrate on early half');
   txtvars=char(txtvars,'  6 RE for validate on late half');
   txtvars=char(txtvars,'RMSE root mean square error of full-period recons and cross-valid');
   txtvars=char(txtvars,'  1 calibration');
   txtvars=char(txtvars,'  2 cross-validation');
   txtvars=char(txtvars,'C1comm  calibration R-sqd after refitting models to a common period');
   txtvars=char(txtvars,'C4 regression coefficients for full-period model, as follows:');
   txtvars=char(txtvars,'   constant, lag 0, lag t-1, lag t-2... lag t+1, lag t+2,....');
   txtvars=char(txtvars,'txtsum text summary of settings for analysis');
   txtvars=char(txtvars,'txtvars  definitions of saved variables');
   txtvars=char(txtvars,'table1, table2, ... -- summary string tables ');
   
   RMSE=RMSE(~Lomit,:);
   REcomp=REcomp(~Lomit,:);
   C1comm=C1comm(~Lomit);
   C4=C4(~Lomit,:);
   NP=NP(~Lomit);
   
   saveset1=' Lomit nslean X2 REcomp RMSE NP C4 C1comm ';
   saveset2=' table1 table2 table3 table4 table5 txtsum txtvars';
   saveset=[saveset1 saveset2];
   % Store results
   eval(['save ' pf5  saveset ';']);
end

fclose all


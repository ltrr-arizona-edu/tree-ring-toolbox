function reglag2
% reglag2:  distributed lag reconstruction of tree-centered reconstructions of pcp 
% CALL: reglag2
%
%****************  IN **********
%
% No input arguments
%
% User prompted to point to input control file with following line by line:
%
% number of tree-ring series (either tree indices or chronologies)
% path\file of .mat file holding the tsm of tree-ring series-- data in X, year
%		vector in yr
% path to .mat files holding the tree-interpolated pcp series
% seascode -- season code of water balance file <oc_se>
% yrgo, yrsp =start, end years of period to consider for modeling
% nlag -- number of + and - lags 
% kmode (1 or 2): exploratory regression mode (1) or reconstruction mode 2)
% path\file of .mat file holding cv  Lomit, which says which, if any,  of the
%    tree sites to omit from use in the subsequent PCA regression. 
%
%******************* NOTES *************
%
% reglag2.m written specifically for San Pedro River Basin pcp reconstruction
%
% Works on multiple tree-ring series
%
% Validates models by split-sample modeling.  Splits available data in
% yrgo to yrsp into two periods, second half largest if total number of 
% years odd.  Calibrates on second, verifies on first. Calibrates on first,
% validates on second.  Then fits entire data period for final model.
% Saves calibration and validation stats for split-sample modeling.  Saves
% calibration stats for full-period model.  In optional fine-detail mode, 
% gives graphs of residuals analysis for full model
%
% Missing data handling.  The period yrgo-yrsp is the longest possible modeling
% period.  If tree or climate series does not cover this whole period, the
% (shorter) common period of climate and tree data is used for modeling
%
% Accepts surrogate data (e.g., made by surrgt1.m, etc).  A gap of NaN might
% therefore exist between the calibration period and earlier data
%
%
% STEPS
%
% * Read input controls; size variables; allocate
% * Set data applicable to all tree-ring series
% * Loop over tree-ring series
%   * Record time-period of tree ring series
%   * Load the pcp file
%      * Identify the usable periods (3) for split-sample modeling
%      * Calibrate on second half, validate on first
%      * Calibrate on first half, validate on second
%      * Full-period modeling
%      * Save key results in file
% * end loop over tree-ring series 
%
%
%
%************************* TABLE SUMMARY ***********************
%
% Need to design table with following cols
%
% 	Site # (* for marked for omission from regl recon)
%		Tree identifier
%		species
%		location
%		elev
%		start, end of good data period
%
% And another with this
%
%		Site #  (* for marked for omission from regl recon)
%		Short name
%		wb variable recon
%		sub-region
%		ARmodel
%		AR var expld
%		SS reg coefs: const -2, -1 0 1 2
%		R-sqd
%		F(p-value)
%		
%
%
%********************* GRAPH SUMMARY
%
% Need to build in capability to plot observed and recon  variable, with
% annotated:  
%
%		Site # and name, species
%		R-sqd, F(p-value)
%		y-variable labeled along y axis
%
%****************  END OF OPENING COMMENTS

a=NaN;


%---------- READ INPUT CONTROLS
[file1,path1]=uigetfile('*.ctl','Input control for reglag2.m');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

% Number of tree ring sites
nsites = str2num(fgetl(fid1)); % number of tree-ring sites

% File of the tree-ring data
pf2 = fgetl(fid1);
eval(['load ' pf2]);
if ~exist('X',1) | ~exist('yr',1) | ~exist('ICD',1);
   error(['X, yr, or ICD does not exist in ' pf2]);
end
[mX,nX]=size(X);
if nsites ~= nX | nsites~=length(ICD);
   error('Specified nsites inconsistent with length of ICD or with col size of X');
end
% store year column for X, and remove year col from X
yrX  = yr;
if length(yr)~=mX;
	error(['length of yr does not match row size of X in ' pf2]);
end
yrsX=[min(yrX) max(yrX)];



% Path to the tree-interpolated pcp files
path3=fgetl(fid1);

% Season code specifying season of pcp grouping
seascode=fgetl(fid1);

%-start end year of period to be considered viable for modeling
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

% Number of positive and neg lags
nlag=str2num(fgetl(fid1));

% Mode (exploratory regression or final filtering)
kmode = str2num(fgetl(fid1));
 

%---------  SET DATA APPLICABLE TO ALL TREE-RING SITES

% -------------Allocate
C1a = repmat(NaN,nsites,1);  % Calib-period R-squared, calib on last half
C1b = C1a; % Calib R-sqd, calib on first half
C1 = C1a; % Calib R-sqd, full model
C2 = C1a; % Calib F-level, overall
C3 = C1a; % Calib prob for F overal
V1a =C1a; % validation rmse, first-half validation
V1b = C1a; % validation rmse, second-half validation
V2a = C1a; % RE
V2b = C1a; % RE

C4 = repmat(NaN,nsites,2*nlag+2); % Coefficients on the various lags
C5 = C4; % 2 standard devs of the coefs

% Build lagged predictor matrix
[XL,yrXL] = lagyr3(X,yrsX,[0 nlag nlag]);
yrXL = (yrXL(1,1):yrXL(1,2))';
% note: yrXL not used in this pgm

% If in exploratory regression mode, truncate rows of XL to hold only the
% years specified by yrgo and yrsp
if kmode==1;
   Ltemp = yrX>=yrgo  & yrX<=yrsp;
   XL = XL(Ltemp,:);
   yrXL = (yrgo:yrsp)';
end

% Set a ones vector and other things needed for long-term reconstruction in
% runs with kmode==2.  If kmode==1, these things are not used, but are 
% nevertheless generated here
[m2,n2]=size(XL);
YL = a(ones(m2,1),ones(nsites,1));  % will hold re-scaled filtered series
onecv = ones(m2,1);

% Clear unneded variables
clear X Ltemp


%----------- LOOP OVER TREE-RING SITES
clc
disp('STARTING LOOP OVER TREE-RING SERIES')
for ns = 1:nsites;
   disp(['Series # ' int2str(ns)]);

	% Get the corresponding chronology site
	crnsite = ICD(ns);

   % Compute rv indexing cols of XL for this series
   mintemp = ns;
   maxtemp = ns + (nlag*2)*nsites;
   icol1 = mintemp:nsites:maxtemp;
   clear mintemp maxtemp;
   
   if ns==19;
      disp('site 19');
   end
   
   
   % Find the years in the yrgo:yrsp period for which no data are
   % missing for the tree-rings at current or lagged values
   Ltemp  = yrXL>=yrgo & yrXL<=yrsp;
   Xtemp= XL(:,icol1);
   % Computer logical pointer to years in Xtemp that all series have data for
   if length(icol1)==1; 
      Lvalid = ~isnan(Xtemp);
   else
      Lvalid =  (~any(isnan(Xtemp')))';
   end
   % Compute vector of years tree-ring data allows modeling for this site
   yrgoodx = yrXL(Ltemp & Lvalid);  
   yrsgoodx=[min(yrgoodx) max(yrgoodx)];
   clear Lvalid Ltemp Xtemp
   
   
   % Get the tree-interpolated monthly pcp file
   file3=['grdt' int2str(crnsite)];
   pf3=[path3 file3];
   eval(['load ' pf3]); % monthly data is in Z, with year in col 1
	if ~exist('Z',1);
		error([pf3 ' does not contain Z']);
	end
	[mZ,nZ]=size(Z);
	zgo = Z(1,1);  zsp = Z(mZ,1);
	
	% Convert monthly pcp to seasonalized pcp
	if strcmp(seascode,'nv-ap');
		begmo = 11; endmo = 4;
	else
		error('Must modify code to reconstruct other than Nov-Apr pcp');
	end
	F = seaspt(Z,begmo,endmo,[zgo zsp],1); % note, the 1 means pcp

   % Put seasonalize data in vector W, year vector in yrW
   yrW = F(:,1);
   yrsW = [min(yrW) max(yrW)];
   W=F(:,2); % drop year col
   
   
   if kmode ==1; % exploratory mode -- not final reconstruction mode
      w = W; % the seasonal pcp series
      
      % Get years for which w has valid data
      Ltemp =~isnan(w);
      yrgoodw = yrW(Ltemp);
      yrsgoodw = [min(yrgoodw) max(yrgoodw)];
      
      % Mark end yrs of period with both clim data and valid lagged data of tree
      % This is the full period for modeling
      yrs3= [max(yrsgoodw(1),yrsgoodx(1))  min(yrsgoodw(2),yrsgoodx(2))];
      yr3 = (yrs3(1):yrs3(2))';
      nyrs3 = length(yr3);
      
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
         
      
      % Pull off the full-period data
      Xkey = XL(LXL3,:);
      wkey = w(Lw3);
      
		% Do the regression
      [I2,I4,stats,c,e,what]=stepr2(Xkey,wkey,icol1);
      
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
      
      if kmode==1;
      %if  Fp>.05;
         %figure(1);
         %plot(yr3,wkey,yr3,what);
         %title(['Series # ' int2str(ns) ': p-value(F) = ' num2str(Fp)]);
         %pause
         %disp(c)
         %disp(I4)
         %pause
       %end
      end
      
   
elseif kmode==2;  % Final filtering mode
      w = W; % the pre-selected water-balance time series
      
      % Get years for which w has valid data
      Ltemp =~isnan(w);
      yrgoodw = yrW(Ltemp);
      yrsgoodw = [min(yrgoodw) max(yrgoodw)];
      
      % Mark end yrs of period with both clim data and valid lagged data of tree
      % This is the full period for modeling
      yrs3= [max(yrsgoodw(1),yrsgoodx(1))  min(yrsgoodw(2),yrsgoodx(2))];
      yr3 = (yrs3(1):yrs3(2))';
      nyrs3 = length(yr3);
      
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
         
      
      % Pull off the full-period data
      Xkey = XL(LXL3,:);
      wkey = w(Lw3);
      
            
      
		% Regression
      [I2,I4,stats,c,e,what]=stepr2(Xkey,wkey,icol1);
      
            
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
      
      
      %************ LONG TERM RECONSTRUCTION
      XTEMP = [onecv XL(:,I2)]; % predictor matrix
      y =  XTEMP*c;  % filtered series
      YL(:,ns)=y;
    end
   
end; % of loop ns over tree ring sites


if kmode==2;
    
   % Re-scale filtered tree-ring series so that variance is proportional to
   % the R-sqd of their single-site regression model
   
   % Use 1700-1961 to get mean and stddev
   L1temp = yrXL>=1700 & yrXL<=1961;
   Z = YL(L1temp,:);
   [mz,nz]=size(Z);
   zmean = nanmean(Z);
   zsd = nanstd(Z);
   
   % Convert YL to zscores
   Y2 = (YL -(repmat(zmean,m2,1))) ./  (repmat(zsd,m2,1));
   
   % Scale by R-sqd
   fact1=C1';
   fact1 = repmat(fact1,m2,1);
   Y2 = Y2 .*  fact1;
   X2=[yrX Y2];
   
   [file5,path5]=uiputfile('resc*.mat','File to hold filtered, scaled series');
   pf5=[path5 file5];
   txt = 'X2 -- tsm of filtered, re-scaled tree-ring series';
   txt=char(txt,'C1 -- R-squared for full-period model');
   txt=char(txt,'C2 -- overall F for full-period model');
   txt=char(txt,'C3 -- p-value for overall F');
   txt=char(txt,'C4 -- regression coefficients, constant first, then neg lags, then +');
   eval(['save ' pf5  ' X2 C1 C2 C3 C4;']);
   
end

fclose all






function regold
% this is temporarily named regold.m while I make new version of reglag1
% reglag1:  distributed lag regression to choose clim variable best related to trees
% CALL: reglag1
%
%****************  IN **********
%
% No input arguments
%
% User prompted to point to input control file with following line by line:
%
% number of tree-ring sites
% path\file of .mat file holding the tsm of tree-ring variables
% path to .mat files holding the water balance seasonal output
% seascode -- season code of water balance file <oc_se>
% yrgo, yrsp =start, end years of period to consider for modeling
% nlag -- number of + and - lags 
% kmode (1 or 2): exploratory regression mode (1) or reconstruction mode 2)
% path\file of .mat file holding cv  iclim and Lomit, which says which of the
%   16 water balance variables to reconstuct.  iclim could be a scalar, 
%   telling that all sites to use same water balance variable, or could
%   be a cv of length equal to number of sites, allowing a different
%   water balance variable to be reconstructed for each site. iclim is not
%   even used if kmode(1)==1. Lomit is logical cv telling which tree sites
%   to omit from use in the subsequent PCA regression.  Lomit and iclim
%   are computed using script file markem1.m
%
%******************* NOTES *************
%
% reglag1.m written specifically to evaluate which of several types of water
% balance variables are best related to tree-ring series. Written for use in
% Sacramento River reconstruction.
%
% Works on 1 tree-ring series, 16 different water-balance variables
%
% Validates models by split-sample modeling.  Splits available data in
% yrgo to yrsp into two periods, second half larges if total number of 
% years odd.  Calibrates on second, verifies on first. Calibrates on first,
% validates on second.  Then fits entire data period for final model.
% Saves calibration and validation stats for split-sample modeling.  Saves
% calibration stats for full-period model.  In optional fine-detail mode, 
% gives graphs of residuals analysis for full model
%
% Missing data handling.  The period yrgo-yrsp is the longest possible modeling
% period.  If tree or climate series does not cover this whole period, the
% subset of periods in common between the two types of data is used for modeling.
%
%
% STEPS
%
% * Read input controls; size variables; allocate
% * Set data applicable to all tree-ring series
% * Loop over tree-ring series
%   * Record time-period of tree ring series
%   * Load the climate file
%   * Loop over the 16 climate variables
%      * Identify the usable periods (3) for modeling
%      * Calibrate on second half, validate on first
%      * Calibrate on first half, validate on second
%      * Full-period modeling
%      * Append key results to files that record results for all 16 vbls
%   * end loop over the 16 climae variables
% * end loop over tree-ring series 
%
%
%
%************************* TABLE SUMMARY ***********************
%
% Need to design table with following cols
%
% 	Site # (* for marked for omission from regl recon)
%		Site name
%		specie
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
% And this for exploratory mode only
%  
%  
%
%********************* GRAPH SUMMARY
%
% Need to build in capability to plot observed and recon wb variable, with
% annotated:  
%
%		Site # and name, species
%		R-sqd, F(p-value)
%		y-variable labeled along y axis





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
if ~exist('X',1);
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


% Path to the watbal output files
path3=fgetl(fid1);

% Season code specifying name of matrix of seasonalize water balance output
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


% File with cv telling which clim variable to reconstruct when kmode==2
pf4 = fgetl(fid1);
eval(['load ' pf4]);
if ~exist('iclim',1);
   error([pf4 ' does not contain iclim']);
end
   
   


%---------  SET DATA APPLICABLE TO ALL TREE-RING SITES

% -------------Allocate
C1a = a(ones(nsites,1),ones(16,1));  % Calib-period R-squared, calib on last half
C1b = C1a; % Calib R-sqd, calib on first half
C1 = C1a; % Calib R-sqd, full model
C2 = C1a; % Calib F-level, overall
C3 = C1a; % Calib prob for F overal
V1a =C1a; % validation rmse, first-half validation
V1b = C1a; % validation rmse, second-half validation
V2a = C1a; % RE
V2b = C1a; % RE

C4 = a(ones(nsites),ones(16,1),ones(2*nlag+2,1)); % Coefficients on the various lags
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

% If in final filtering mode, set a ones vector and other things
[m2,n2]=size(XL);
YL = a(ones(m2,1),ones(nsites,1));  % will hold re-scaled filtered series
onecv = ones(m2,1);


% Clear unneded variables
clear X Ltemp


%----------- LOOP OVER TREE-RING SITES

if kmode==1;
   [file6,path6]=uiputfile('resar*.mat','Store exploratory results here');
   pf6=[path6 file6];
   %fid6=fopen(pf6,'w');
   %fprintf(fid6,'%s',head6);
end

clc
disp('STARTING LOOP OVER TREE-RING SITES')
for ns = 1:nsites;
   disp(['Site # ' int2str(ns)]);
   % Compute rv indexing cols of XL for this site
   mintemp = ns;
   maxtemp = ns + (nlag*2)*nsites;
   icol1 = mintemp:nsites:maxtemp;
   clear mintemp maxtemp;
   
   % Compute first, last year for which the submatrix of current and lagged
   % tree rings for this site have no missing data
   Xtemp = XL(:,icol1);
   Ltemp = (~any(isnan(Xtemp')))';
   yrgoodx = yrXL(Ltemp);
   yrsgoodx=[min(yrgoodx) max(yrgoodx)];
   
   % Get the water-balance file
   file3=['wbout' int2str(ns)];
   pf3=[path3 file3];
   eval(['load ' pf3]);
   clear datout; % do not need this huge cell variable
   % Put the matrix of climate variables in W
   eval(['W = ' seascode ';']);
   
   % Get years for water balance matrix
   yrW = W(:,1);
   yrsW = [min(yrW) max(yrW)];
   W(:,1)=[]; % drop year col
    
    %  vbl 1 -- pcp
    %   vbl 2 -- tmp
    %   vbl 3 -- RO  runoff
    %   vbl 4 -- Z index
    %   vbl 5 -- X  PDSI
    %   vbl 6 -- W  soil moisture (mid month)
    %   vbls 7-16 SI  stress index, which is P-cPE, where c = .1 ,.2, ... 1.0
    %
    
    if kmode ==1;
              
       
       
       
    for nv = 1:16; % Loop over water balance variables
      w = W(:,nv); % a water-balance time series
      
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
      
      [I2,I4,stats,c,e,what]=stepr2(Xkey,wkey,icol1);
      
      
      
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
      
      if kmode==1;
      %if (nv==5 | nv==1) & Fp<=.025
         %figure(1);
         %plot(yr3,wkey,yr3,what);
         %title(['Series # ' int2str(ns) ': p-value(F) = ' num2str(Fp)]);
         %pause
         %disp(c)
         %disp(I4)
         %pause
         % end
      end
      
           
            

   end; % of loop nv over water balance variables
   
      
elseif kmode==2;  % Final filtering mode
      w = W(:,iclim(ns)); % the pre-selected water-balance time series
      
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
      
      [I2,I4,stats,c,e,what]=stepr2(Xkey,wkey,icol1);
      
      
      
      % Store results
      C1(ns,1)=stats(1); % R-sqd for full-period model
      C2(ns,1)=stats(3); % overall F for equation
      
      % P value for overall F   
      Fcompute=stats(3);
	   Fp=1-fcdf(Fcompute,(length(c)-1),(nyrs3-length(c)));
      C3(ns,1)=Fp;
      
      % Coefficients
      C4(ns,1,1)=c(1);
      C4(ns,1,I4+1)=c(2:length(c));
      
      if ns==44;
         disp('wow, site 44');
      end
      
         
      %************ LONG TERM RECONSTRUCTION
      XTEMP = [onecv XL(:,I2)]; % predictor matrix
      y =  XTEMP*c;  % filtered series
      YL(:,ns)=y;
      
      
    end
 
   
end; % of loop ns over tree ring sites


if kmode==2;
   % trim off unneeded dimensions from statistics matrices
   C4=C4(:,1,:);
   C1=C1(:,1);
   C2=C2(:,1);
   C3=C3(:,1);
   
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
   
   [file5,path5]=uiputfile('tsm*.mat','File to hold filtered, scaled series');
   pf5=[path5 file5];
   eval(['save ' pf5  ' X2;']);
   
elseif kmode==1;
   txt = 'C1 -- R-squared';
   txt=char(txt,...
      'C2 -- overall F',...
      'C3 -- p-value for overall F');
   eval(['save ' pf6  ' C1 C2 C3 txt ;']);
   
   
end

fclose all


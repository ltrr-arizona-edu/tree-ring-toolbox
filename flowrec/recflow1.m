function recflow1
%recflow1: flow reconstruction by PCR, with cross-validation
%CALL: recflow1
%
% Meko 10-15-97
%
% Revised 7-31-98:  to return matrix of calib-period recons for all sub-period
%   models
%
% Features: 
% -different recon models for different time periods depending on tree-ring coverage
% -cross-validation tests sensivity of results to number of PC's retained as potential
%  predictors and to number of those used as final predictors
% -follows application of reglag1.m to filter and scale tree-ring chronologies
%  according to strength of local climate signal for water-balance variables
%
%************ IN *****************************
%
% No input arguments.  User is prompted to select .mat files with input
%
% File 1 <ts*.mat>  file with filtered, scaled tree-ring site chrons in matrix
%		X.  Typically, this matrix built using reglag1.m.
% File 2 <*.mat>  2-col .mat file with year in col 1 and the variable
%		to be reconstructed (e.g., water -year flow of Sacramento River) in 
%		col 2 
% File 3 <recin*.mat> specifications for the reconstruction model:
%
%		Assume np reconstruction sub-periods and ns tree-ring
%		sites.  Variables in file 3 are are as follows:
%
%    kmode (1 x 1)i running mode (see notes)
%        ==1 reconstruction-1: automatically fits 'reasonable' model  
%        ==2 exploratory cross-validation; points out possibly better model
%        ==3 reconstruction-2: final fitting of specified model
%		yrsC (1 x 2)i  start, end year for calibration 
%		kopt (1 x 1)i  options for modeling:
%			kopt(1) -- debugging or cranking 
%						== 1	diagnostic feedback -- debugging mode
%						== 2	vital feedback only -- cranking out results
%			kopt(2) reserved
%			kopt(3) reserved
%			kopt(4) -- default predictor selection option for PCR
%					1 -- use all potential predictors as predictors
%					2 -- use only those potential predictors (PC scores) whose
%					  regression coefficients are significant at the significance level
%				     tcrit by a t-test (Mardia et al. 1979)
%    mineigs (1 x 1)i minimum number of eigenvectors for dataset reduction
%			by discarding predictor eigs.  If a given time-period model has
%			only a few sites, and so only a few PC's, why not use them all as
%			potential predictors in the PCR models?  Setting for mineigs lets
%			you do this.  For example, setting mineigs to 7 would mean as long
%			as the number of sites in region is 7 or fewer, you would let all 
%			PC's be potential predictors.  Otherwise, would separate wheat from
%			chaff using, say, an eigenvalue-1 threshold
% 	   tcrit (1 x 1)r signif-level threshold for entry of predictors in the PCR
%			models (typically, set tcrit=0.95). Used only if kopt(4)==2
%		evrule (1 x 1)r  default rule for how many predictor pc's to retain as
%			potential predictors for PCR  (see notes) 
%				evrule == 1   retain all PCs
%				evrule == 2   discard PCs whose eigenvalue smaller than ave
%           evrule == 3   retain enough pcs to explain fract1 of variation
%	   nkeep(1 x np)i  # of pc's to retain as potential predictors
%		npred (1 x np)i  or number of retained PC's to include as
%				predictors in the final equation
%     fract1 (1 x 1)r retain enought PCs to explain this decimal fraction
%        predictor variance
%
% file5 -- .mat input file with the string info on tree-ring file names,
%   site names, state, and species <sitenms>, with info in Snames
% file6 -- .mat input file with Lpergo and Lperiod, which define the 
%   subperiods for regional reconstruction. These variables were originally 
%   generated and stored by script file subper1.m
%
%	  Lperiod (ns x np)L  cross-reference logical matrix indicating which
%        sites are to be included in which sub-period models.  (1) means
%        include, (0) means do not include.  Lperiod built previously using
%        say, subset1.m
%    Lpergo (1 x np)i start year of each sub-period
%%
% file7 -- .mat input file with Lomit -- which series not to use
%		Lomit (ns x 1)L  whether site is to be omitted (1) or not (0) from
%			the regional analysis.  This was determinted from reglag1.m results.
%       Lomit allows user to omit sites marginally correlated with climate
%
%************  NOTES *********************
%
%
% This reconstruction method was designed for application to reconstruct a
% regional time series (e.g., flow of a river) from a spatially distributed
% network of tree-ring sites whose trees might vary site-to-site in which
% climate variable governs growth, and how strongly. Function recflow1.m is
% meant to be run following application of reglab1.m, which generates a
% matrix of time-filtered tree-ring series whose variance is scaled 
% proportional to the local climate signal in the tree-ring series.  Function
% recflow1.m was written to satisfy needs of a good reconstruction method for
% annual (wy) flow of the Sacramento River from 57 tree-ring chrons distributed
% over Calif, Oreg and Nev.
% 
% Cols of input tree-ring matrix X are for individual tree-ring sites. 
% Variables have been generated by reglag1.m, and have variance equal to
% the decimal proportion of variance of a selected water-balance variable
% (e.g., water-year PDSI,  water-year precip) interpolated to the tree site
%
% kmode.  Running mode for function.  Exploratory cross-validation is used
% for diagnostic data on how many PC's to keep as potential predictors, and
% on how many of those to use as actual predictors.  Reconstruction-2 mode 
% lets you set those quantities for each sub-period model.  Reconstruction-1
% mode lets you determine those quantities automatically using various rules.
% Smart way to proceed is to (1) run recon-1 mode to get feel for range of
% reasonable settings, and for first look at what accuracy of recon might
% expect, (2) set fract1 and run in exploratory cross-validation mode to
% check where overfitting definitely becomes a problem, and (3) specify
% nkeep and npred and run final model in recon-2 mode.
%
% kopt(1).  Debugging vs cranking.  Debugging mode gives extra graphical
% 	and screen text output.  Otherwise runs same as cranking mode.
%
%
%
% PCA is done on covariance matrix in this program.  Covariance matrix is used 
% because want to retain influence of variance differences from site to site
% in the filtered, scaled tree-ring series in input matrix X
% 
%
%**********************  STEPS **********************************
%
% Check input
% Build year pointers
% Store calib-period means, std devs, of predictand
% Compute applicable recon year segment for each sub-period
% Allocate to store model output
% Loop over subperiods, in kmode 1,2, or 3
%
% 	In kmode 1.  Uses rules to pick reasonable size model. Generates
% 	a reconstruction. Cross-validates using leave-1-out strategy. 
% 	Yields reconstruction and statistics.
%
% 	In kmode 2.  Computes cross-validation statistics (rmse, RE) for
%	a range of settings of number of PCs to keep as potential predictors,
%	and number of those to use as actual predictors in PCR.  The matrices
%	of stored rmse and RE validation stats for each sub-period can be used
%	as a guide along with results of run with kmode==1 to choose final
%	settings for model specification nkeep and npred.  No long-term 
%    reconstruction generated in kmode 2
%
% 	In kmode 3.  You specify number of PC's to keep as potential 
% 	predictors and number to use as final predictors for each sub-period.
% 	Output just like that for kmode 3.  This is usually the mode used
% 	to get final reconstruction, as it allows user intervention in 
% 	building in info gained from runse in kmode 3 and 1. 
% End Loop
%
% If one of the recon modes, 
%		assemble final recons and error bars from sub-period recons 
%		analyze residuals and build residuals analysis codes
%		build summary tables and figures (see below).
% Else exploratory cross-validation mode
%		Graph results.
% Endif
%
%
%******************** SUMMARY TABLES AND FIGURES *************
%
% 
% Exploratory cross-validation mode (kmode==2)
%
%  Color map of RE values for various nkeep,npred combinations
%  Color map of rmse...
%
%
% Reconstruction modes: kmode 1 or kmode 3
%
%		tsp of reconstruction with error bars (for screen viewing andzooming)
%		tsp (hq) of reconstruction on one page, 200 yrs per axis, with error bars
%		Landscape-oriented table summarizing models; cols as follows
%			Model no.
%			Period for model
%			Number of variables:  tree-ring chrons, pot predictors, final preds
%			Rsqd for calib of full model
%			F (pvalue) for full model
%			Reconstruction error (+- corresponds to 95% error bars)
%			RMSE (2 cols)
%				Calibration
%				Cross-validation
%			RE statistic
%			Residuals code:  a=autocorr resids, b=time trend in resids
%					c=resids variance function of actual dependent variable
%					d=resids highly non-normal
%
%		Ascii file:  year, estimated predictand, std error of prediction,
%				sub-period model #



%*********** CHECK INPUT *************************


%--------------- TIME SERIES MATRIX OF PREDICTORS

a=NaN; 

[file1,path1]=uigetfile('*.mat','.mat file with tsm of predictors in X2');
pf1=[path1,file1];
eval(['load ' pf1 ' X2;']);
if ~(exist('X2')==1);
   error([pf1 ' does not contain X2']);
end
[mX,nX]=size(X2);

%-------  MATRIX OF STRING SITE NAMES
[file5,path5]=uigetfile('site*.mat','.mat file site text info');
pf5=[path5,file5];
eval(['load ' pf5]);
if ~(exist('Snames')==1);
   error([pf5 ' does not contain Snames']);
end

%-------  Lomit  -- which tree-site series to not use in the regional model
[file7,path7]=uigetfile('clim*.mat','input .mat file with Lomit');
pf7=[path7,file7];
eval(['load ' pf7 ' Lomit;']);
if ~(exist('Lomit')==1);
   error([pf7 ' does not contain Lomit']);
end



%------------  PREDICTAND 
%
% Assumes a .mat file with matrix y: year in col 1 and predictand in col 2 is available
% 
[file2,path2]=uigetfile('*.mat','.mat file with predictand y');
pf2=[path2 file2];
if exist('y')==1 | exist('Y')==1; 
   error(['y or Y already exists before loading ' pf2]);
end
eval(['load ' pf2]);
if ~(exist('y')==1)  & ~(exist('Y')==1);
   error([pf2 ' does not contain y or Y']);
end
if exist('Y')==1;
   y=Y;
   clear Y;
end

[my,ny]=size(y);
if ny~=2; 
   error('y must be 2-col matrix');
end

% Allow for log transform of predictand
buttonname=questdlg('Log10 Transform the Predictand?');
switch buttonname
case 'Yes';
   y(:,2)=log10(y(:,2));
case 'No';
otherwise
end

%-------------- INPUT Lperiod and Lpergo
[file6,path6]=uigetfile('pty*.mat','input file with Lperiod and Lpergo');
pf6=[path6 file6];
eval(['load ' pf6 ' Lperiod Lpergo;']);



%------------ MISCELLANEOUS CONTROL

[file3,path3]=uigetfile('recin*.mat','recin?.mat file with misc input');
pf3=[path3 file3];
eval(['load ' pf3]);

% Check for consistency of dimensions of X2 tsm with misc matrices
[mtemp,ntemp]=size(Lperiod);
if nX-1 ~= mtemp;
   error('Rows in Lperiod should match cols  (minus year) in X2)');
end
np = ntemp; % number of periods
ns = nX-1; % number of sites in X2


%**************  BUILD YEAR POINTERS **************************


% Compute start, end years of calib and valid periods
% of A-models and  B-models, and the calibratino period for the full model C
% A-model -- calibs on most recent period, validates on earliest period
% B-model --  calibrates on early years, validates on recent years
yrC = [yrsC(1):yrsC(2)]'; % vector for full calib period
nyrsC = length(yrC);

% Get year info off full tree-ring matric X2, then strip year col off X2
yrX2 = X2(:,1);
X2(:,1)=[];

% Get year off predictand 2-col matrix.  Then strip the year off.
yry = y(:,1);
y(:,1)=[];

% Make pointers to full calib period (period 'C')
LCX2 = yrX2>=yrsC(1) & yrX2<=yrsC(2); % in full tree-ring matrix X2
LCy = yry>=yrsC(1) & yry<=yrsC(2);  % in predictand series

%*********  COMPUTE APPLICABLE YEARS FOR RECONSTRUCTION SEGMENTS

yrtemp = Lpergo; 
yrtemp(1)=[];
yrtemp=[yrtemp yrsC(1)];
yrssub = [Lpergo; yrtemp-1];  % row 1: start year. row 2, end year 


%*********** ALLOCATE TO STORE MODEL OUTPUT

Y2 = a(ones(mX,1),ones(5,1)); % to hold year, recon value rmse
Y2(:,1)=yrX2;
Ecv = repmat(NaN,nyrsC,np);  % to hold cross-validation calibration period residuals
%       for each sub-period model

NV = a(ones(np,1),ones(3,1)); % number of variables in models
% 1,2,3 -- number of tree-ring sites avail, # potential predictors, # selected
Rsq = a(ones(np,1),:); % calib R-sqd 
rsq = a(ones(np,1),ones(2,1)); % valid r-sqd for A, B models
F = a(ones(np,1),:); % overall F for final model
Fpval = F; % p-value for overal F
ebar2 = F ; % +- this for 5% confid limits
RMSE = a(ones(np,1),ones(2,1)); % root mean square error for (1) calib and
%		(2) cross validation
RE = a(ones(np,1),:); % reduction of error stat, A model, B model
blnk6=blanks(6);
rcode = repmat(blnk6,np,1); % residuals analysis code

% Revision 7-31-98
YPC = repmat(NaN,nyrsC,np);
% Revision end 3-31-98


%********************** LOOP OVER SUB-PERIODS

for n = 1:np;
   
   % Applicable reconstruction years
   yrsapp = (yrssub(:,n))';
   
   % Pointer to the applicable rows of  reconstruction data in X2
   Lsub = yrX2>=yrsapp(1) & yrX2<=yrsapp(2);
   
   disp(['Starting analysis for period: ' int2str(yrsapp(1)) '-' int2str(yrsapp(2))]);
   
   % Which of the ns sites to use?
   Luse =~Lomit & Lperiod(:,n);
   nsites = sum(Luse);
   if nsites==0;
      error(['No sites qualify for period  starting ' int2str(yrsapp(1))]);
   end
   
	% Get time series of predictor data for applicable reconstruction years
	 Xlong = X2(Lsub,Luse);
   
	% Get time series of calibration data
  X3 = X2(LCX2,Luse); % predictor
  y3 = y(LCy); % predictand
   
  % Set options for call to pcarec1
  kopta(1)= evrule; % for reducing predictor data to potential predictors
  kopta(2)=1; % reserved
  kopta(3)=kopt(4); % for selecting final predictors from potential predictors
   
  datin{1}=X3; % calib period predictor data
  datin{2}=y3; % calib period predictand series
  datin{3}=Xlong; % applicable reconstruction years 
   
  setmodel{1}=nkeep(n); % retain this many predictor PCs
  setmodel{2}=npred(n); % use this many retained PCs as predictors in regression
   
  specs{1}=mineigs;
  specs{2}=fract1;
  specs{3}= tcrit;
  specs{4} = kopta;
   
  % Estimate and cross-validate reconstruction model
  [recout,cvout]=pcarec1(kmode,datin,setmodel,specs,Snames(Luse,:));
   
   if kmode==1 | kmode==3; % unpack cell data from reconstruction
            
     Rsq(n)=recout{3}; %R-squared for full-model calibration
     RMSE(n,:)=recout{4};
     RE(n)=recout{5};
     NV(n,:)=recout{7}; % number of tree-ring sites for PCA, number of PCs
     % passed to regression as potential predictors; number used in final
     % model
     
     Y2(Lsub,2)=recout{2};  % reconstructed output for pre-calibration segment
     % Slap on the calibration and cross validation rmse
     Y2(Lsub,3)=RMSE(n,1);
     Y2(Lsub,4)=RMSE(n,2);
     
     % Slap on calibration R-squared
     Y2(Lsub,5)=Rsq(n);
     
     % Store cross-validation residuals for calib period
     Ecv(:,n)=recout{10};
     
     % Revision 7-31-98
     YPC(:,n) = recout{1};
     % Revision 7-31-98 end

          
  else; % unpack cell data from exploratory cross-validation
     NK{n}=cvout{1};
     S1{n}=cvout{2}; % RMSE for various combinations of PCs retained and how
     %  many of thos used as predictors
     S2{n}=cvout{3}; % RE for various combinations of PCs retained and how
     %  many of thos used as predictors
          
  end
  
   
end; % of for n=1:np

if kmode==1 | kmode==3;
   txtsave='Output from recflow1.m run ';
   txtadd = ['Results from a kmode = ' int2str(kmode) ' run'];
   txtsave=str2mat(txtsave,txtadd);
   txtadd='Saved Variables Defined as follows:';
   txta1 = 'Y2 (5 cols): year,reconstructed value,RMSEcal,RMSEval,R-sq calib';
   txta2 = 'Rsq: R-squared for full-model calibration';
   txta3 = 'RMSE (2 cols): calib and validation RMSE';
   txta4 = 'RE: RE cross validation';
   txta5 = 'NV (1 x3)i: numbers of variables, for summary table';
   txtsave=str2mat(txtsave,txtadd,txta1,txta2,txta3,txta4,txta5);
   txta1 = 'y (cv): actual predictand for period yry';
   txta2 = 'yry (cv): year vector for y';
   txta3 = 'Ecv  (np cols)  cross-validation resids for calib period, all models';
   txta4 = 'YPC: calibration-period reconstructed data for each subperiod model';
   txtsave=str2mat(txtsave,txta1,txta2,txta3,txta4);
   txtwhat = txtsave;
   clear txtadd txtsave txtadd txta1 txta2 txta3 txta4;
   
   [file4,path4]=uiputfile('rfo*.mat','Store output here');
   pf4 = [path4 file4];
   set1 = ' txtwhat Rsq RE RMSE Y2 NV y yry Ecv YPC';
   eval(['save ' pf4 set1]);
   
elseif kmode==2;
   txtsave='Output from recflow1.m run';
   txtadd = ['Results from a kmode = ' int2str(kmode) ' run'];
   txtsave=str2mat(txtsave,txtadd);
   txtadd='Saved Variables Defined as follows:';
   
   txta1 = 'NK{i}: number of predictors corresp to fract 1 for ith subset model';
   txta2='S1{i}: matrix of RMSE for various combinations of PCs retained and used';
   txta3='S2{i}: matrix of RE for same';
   txtwhat = str2mat(txtsave,txtadd,txta1,txta2,txta3);
   
   clear txtadd txtsave txtadd txta1 txta2 txta3;
   
   [file4,path4]=uiputfile('rfo*.mat','Store output here');
   pf4 = [path4 file4];
   set1 = ' txtwhat NK S1 S2 ';
   eval(['save ' pf4 set1]);
end





fclose all;

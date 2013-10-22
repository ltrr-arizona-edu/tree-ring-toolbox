% spatrec1.m
%
%  D Meko, 12-18-95
%
% Climate reconstruction from time-varying subsets of tree
% indices. Lags included as predictors in single-site
% reconstruction models are t-1,t,t+1.  Time series of one or
% more "regional" climate variables are reconstructed.
% 
% Station time series of local climate variables are averaged
% to form tree-centered local climate series.
%
% Regional time series (the objective of reconstruction) are 
% either read in or formed by averaging the local climate variable
% over stations.
%
% A local climate series is first reconstructed for each tree
% using lagged regression on the tree index.  The reconstructions
% for the various trees are called "single-site reconstructions"
% (SSR).  
%
% The SSR's are optionally averaged over specified trees to 
% reduce dimensions of the data set.
%
% Subsets of trees living in different periods of the tree-ring
% record are identified.
%
% Principal components regression (PCR) models are estimated
% for each period.  The predictors are PC scores of either the
% SSR's for active trees, or for SSR's averaged over trees.
%
% The final reconstruction is generated from the PCR models, and
% the reconstruction is cross-validated.
%
%******************* INPUT *****************************
%
% T (n x nt)r  tree indices, n years, nt trees; missing
%		values NaN; no year col
% C (mc x nc)r "local" climate matrix, mc years, nc stations;
%		missing values NaN 
% ICD (m1 x n1)i index pointer to columns of local climate matrix
%		C indicating which stations should be averaged to
%		form the tree-centered climate series to be used as
%		predictands for SSR models. Rows of ICD correspond to trees,
%		columns to climate stations. The Column-size  of ICD  depends
%		on the maximum number of stations averaged for any tree.
%		ICD is right-filled with NaN.
%		Example: if climate stations 3,6 and 9 apply to tree #1,
%			and the maximum number of climate stations for any tree 
%			is 6, row 1 of ICD would be:  
%			3  6  9  NaN  NaN NaN
%			The single-site reconstruction for tree#1 would be PPT
%				or whatever averaged over stations 3,6 and 9
% ICY (m3 x n3)i index pointer to columns of station climate matrix
%		C indicating stations to be averaged to form regional
%		climate series. Row size m3 depends on 
%		how many regional series are wanted.  In simplest case of
%		averaging all station climate series into one region, m3=1.
%		Column size varies as in ICD. 
%		*** Note that ICY is not needed if the regional climate series
%		are input directly (kopt(4)==2)
% YY (mYY x nYY)r 
%		Regional climate matrix.  Only if kopt(4)==2, is this matrix
%		necessary. Otherwise, the regional climate series are formed
%		as averages over the local climate series. An example where 
%		YY would be read in is reconstruction of gridpoint 700 mb
%		height.  The local climate variable for this example might
%		be station precipitation for the winter season.
% IDX (m4 x n4)i index pointer to columns of the tree matrix T
%		indicating trees the SSR's should be averaged over to reduce the
%		column size of the predictor data before PCR modeling.
%		The row dimension m4 is <= the number of trees nt.  Likewise
%		the column dimension n4.  In the simplest case of no averaging
%		over trees, IDX is a column vector [1 2 3 ... nt]'. IDX is
%		a convenient way to delete individual tree's SSR's from
%		use in the PCR models.  
% YRS (4 x 2)i or  (5 x 2)i  beginning and ending years of:
%		row 1: T, the tree-index matrix
%		row 2: C, the local-climate matrix
%		row 3: calibration period for the model; be sure to 
%			truncate end year from end year of T by number of +
%			lags in single-site reconstruction
%		row 4: reconstruction period: should not overlap calib period; and
%			allow for lag needs
%		row 5: YY, the specified regional climate matrix;  row 5 needed
%			only if kopt(4)==2
% Fcrit (1 x 1)r desired signif level for overal F statistic of SSR
%		models.  Typically, set Fcrit=0.95.  Fcrit is 1 minus
%		the "p-value" for the test. See Weisberg (1985, p. 49).  Note
%		that a constant is in the models.  Thus, say a model has 53 yr
%		for calib, and predictors entering at t and t-1. Then the 
%		relevant F distribution is F(2,50), following Weisberg. Note
%		that Weisberg uses the p-value rather than the signif level in
%		his discussion.  Spatrec1.m sends
%		a warning message to the screen if the computed F-level is not
%		significant at the level Fcrit.  Such information might be used 
%		in exploratory analysis to delete insensitive trees from
%		consideration in the reconstruction.  Note that Fcrit is not
%		used to determine whether additional lags should enter the
%		SSR models.  Variables entry is governed by adjusted
%		R-squared, which must increase if the model is to be expanded
% lg (3 x 1)  shifting and lagging parameters
%		1 number of years to shift T (the delay; always negative)
%		2 number of negative lags to include on T; 
%		3 number of positive lags to include on T;
%			For now, spatrec1.m works only with
%			lg(2)==1 and lag(3)==1;  in other words, a +-1 lag model
% mineigs (1 x 1)i minimum number of predictor amplitudes to use
%		in PCR models. 
%		If model has this or fewer predictor variables available,
%		all predictor PC's are included as potential predictors,
%		not just those retained by following an eigenvalue-1
%		threshold
% tcrit (1 x 1)r signif-level threshold for entry of predictors
%		in the PC regression models (typically, set tcrit=0.95)
% kopt(1 x 4)i options
%		kopt(1)-- rule for retention of important predictor eigs
%			1 -- use all
%			2 -- eigenvalue of 1 cutoff
%		kopt(2) -- like kopt(1), except for predictands
%		kopt(3) -- predictor selection option
%			1 -- use all potential predictors as predictors
%			2 -- use only those potential predictors (PC scores) whose
%				regression coefficients are significant at 
%				the significance level tcrit% by a t-test 
%				(Mardia et al. 1979)
%		kopt(4)-- method for regional climate series
%			1 -- computed by averaging over local-climate stations
%			2 -- specified in YY
%
%************************  NOTES ********************************
%
% Two stages of regression are used.  The first stage is MLR
% of a local climate time series on lagged values of a tree index.
% This step is done separately for each tree.  The regression equation
% is used to generate single-site reconstructions (SSR's).  The SSR
% regression statistics indicate strength and dynamics (lagged
% response properties) of tree indices, and can be used in 
% exploratory analysis to screen out insensitive trees.  Each SSR
% is a time-filtered version of the tree index, with filter weights
% estimated to maximize the explained local-climate variance.
% Another way of looking at the SSR operation is as a deconvolution
% of the climate signal from tree indices.
%
% The second stage of regression is principal components regression 
% PCA of the regional climate variables on the previously generated
% SSR's.  This stage is spatial weighting of the individual tree
% tree climate signal. PCA is intended to reduce data dimensions,
% which are likely to be cumbersome for typical tree-ring problems.
% No lags are included in the PCR model because the lagged response
% has ideally been deconvoluted in the earlier SSR step.
%
% PCA weighting of tree indices essentially bypasses
% computation of site chronologies.Some tree grouping other than 
% that used in forming site chronologies might be more appropriate 
% for capturing the climate signal in the region.  For
% example, several runoff-sensitive trees might be distributed
% over several chronology collection sites in a river basin. 
%
% An important part of the modeling procedure is the calibration and
% cross-validation of separate reconstruction models for various 
% time periods in the tree-ring record depending on the trees alive
% at the time. This approach allows the reconstruction to extend back
% to the start year of the oldest tree at any site, and allows the
% time-variation of reconstruction accuracy to be estimated more 
% directly than by previously existing methods.
% The tree indices in T should have been prewhitened beforehand 
% (e.g., by ARMA modeling) to remove persistence.  Prewhitening will
% reduce problems of multicollinearity in the SSR predictors.
% A possible approach is to use the "residual" tree indices as output
% by ARSTAN or an equivalent program
%
% A restriction on the tree-ring and climate data is that the
% calibration-period portions of the data matrices must contain no
% missing values (NaN's). NaN's are allowed outside the specified 
% calibration period. 
%
% Some allowance for intermediate reduction of the tree-ring data
% set is made through the IDX input matrix.  Without such reduction,
% the number of predictor variables to be reduced to PCR predictors 
% is the same as the number of trees in T.  With reduction, SSR's
% are averaged over trees before identifying subsets of predictors
% for various time periods for the PCR models.  The trade-off is a
% possible truncation of years to the common period for any group of
% trees pointed to by IDX.  For example, if tree 5 covers 1600-1995
% and tree 6 covers 1710-1995, specifying averaging over trees 5
% and 6 will lead to a usable SSR covering only 1710-1995.
%
%
% Limitations. (1) the model could become unwieldy for large-scale
% reconstructions with possibly hundreds of tree-ring chronology
% sites and thousands of trees (yet to be programmed is a site-
% chronology strategy % for data reduction)
% (2) all tree series must completely overlap the calibration period
% with no missing data
%
%
%************************ OUTLINE OF PROGRAM STEPS **************
%
% STEP-1 Check input arguments, size and allocate variables
% STEP-2 Make site-average and regional average climate series
% STEP-3 Generate single site reconstructions
% STEP-4 Optional averaging of SSR's over trees
% STEP-5 Identify different SSR subsets for various time periods
% STEP-6 Generate spatial reconstructions
% STEP-7 Cross-validate reconstruction by leave-5-out method
%
%
%***************** KEY SAVED VARIABLES  **********************
%
% The following variables are saved in a .mat file.  The user is 
% prompted for the file name.  All of these variables will 
% be in the .mat file only if the user specifies (also 
% screen-prompted) to proceed with cross-validation.  Otherwise
% a subset of variables is saved.
%
%Cwgt      ICD       PEbar     RSQ2      hf        nmods     
%DHS       ICY       Pyrs      RSQ3      kcv       tcrit     
%DHkey     IDX       R2        STATS     kopt      yr3       
%DS        L6        R2ss      YO        lg        yr4       
%EV        Mcv1      RE        YP        mineigs   
%EVcv      NV        REss      e         modnum    
%Fcrit     P2        RSQ1      filesv    modyrs    
%
% Definitions below make use of the following dimension variables
%   mTC -- number of years in calibration period, assumed to be
%		the same for all SSR models and PCR models
%   mrec -- number of years in reconstruction period; the 
%		reconstruction period as defined here does not include the 
%		calibration period
%	 nt --  number of trees in input matrix T
%   q -- number of gridpoints or regions to be reconstructed;
%		q also equals the total number of predictand PC's, 
%		before optional reduction by an eigenvalue-1 criterion
%	 qq -- number of retained predictand PCR predictors (qq<=q)
%   nmods -- number of PCR models in the spatial regression
%
% Cwgt (4 x nt)r  regression coefficients for the SSR models
%		rows 1-4: regression constant, coefficients at 
%		lags zero, -1, +1 years from climate year
% DHS (mTC x nt)r SSR's, calibration period
% DHkey (mTC x nt)r SSR cross-validation estimates, calibratin period
% DS (mTC x nt)r observed tree-weighted local climate time series
%		for each tree.  These time series are the predictands for the
%		SSR models. (Depending on IDX, predictands might be averaged
%		SSR's over columns of DS)
%		by IDX
% EV (nmods x 1)r "explained variance" statistic for PCR models 
%		EV takes into account the PCR explained variance of the 
%		predictant PC scores and the fractional variance of the 
%		original predictand variable accounted for by the PC scores.
%		Defined as in eqn 17, Cook et al. 1994 
% EVcv (mTC x nmods)r Like EV, but for calibration of each leave-5-
%		out version of the PCR models.  Note that nmods is column
%		size for EVcv, but row size for EV
% Fcrit (1 x 1)r  critical significance level for overall-F level of
%		SSR models *** SEE INPUT
% ICD (m1 x n1)i index pointer to columns of station climate matrix
%		C indicating which climate stations should be averaged in
%		to form the local-climate predictands for
%		 each tree.  *** SEE INPUT
% ICY (m3 x n3)i index pointer to columns of station climate matrix
%		C indicating stations to be averaged to form regional
%		climate series. *** SEE INPUT
% IDX (m4 x n4)i index pointer to columns of the tree matrix T
%		indicating trees the SSR's should be averaged over to reduce the
%		column size of the predictor data *** SEE INPUT
% L6 (mL6 x nL6)L  indicates trees or sets of trees used in
%		used in each PCR model. A row of L6 corresponds to a PCR
% 		model.  A column corresponds to a row of IDX, which in
%		turn points to a particular tree or group of trees.
%		The entries, "1" or "0" mean "use" or
%		"do not use" the SSR based on that tree (or the tree-
%		averaged SSR, depending on IDX).
%		The number of rows in L6 equals the number of PCR 
%		models. The number of columns equals row-size of IDX. 
% Mcv1 (mTC x nt) calibration-period means of the SSR's
%		for the leave-5-out cross-validation periods.  Note
%		that differences along a column show how unstable the
%		mean is as a function of subset
% NV (nmods x 5)i summary information on predictor and predictand
%		data sets the PCR models (rows of NV). Columns of NV are
%		[q qq p pp nused] number of predictands in pool, number of
%		predictands retained, number of SSR's available, number of
%		PC's of the SSR's retained as potential predictors, and 
%		number of PC's used as predictors in the final model
% P2 (nmods x 2)i  start and end year of period for each spatial model
% PEbar (nmods x q)r root-mean-square error of prediction, from 
%		cross-validation.  Each row applies to a spatial model. Each
%		column to a gridpoint or climate region. Computed as follows:
%		1. compute the predicted error sum of squares (PRESS) from the
%			leave-5-out cross-validation
%		2. for the i th  model, PEbar(i,:) = sqrt(PRESS/mTC); 
% 		Weisberg, top, p. 230 calls PEbar a "sensible estimate of 
%		average prediction error"
% Pyrs (nt x 2)i start and end years for single-site reconstructions;
%		each row corresponds to a tree, as ordered in the columns of
%		the input matrix of tree indices T
% R2 (nt x mTC)r  calibration R-squared for the 
% 		leave-5-out SSR models at each of the nt trees.
% R2ss (nt x 1)r squared correlation coefficient between observed
%		and predicted local climate variable for leave-5-out 
%		cross-validation of SSR models. R2ss is a verification
%		statistic.
% RE (nmods x q)r  reduction-of-error statistic for the PCR models.
%		row-size equals number of models;  column size equals number
%		of gridpoints or regions reconstructed.  RE is computed 
%		directly from the RMSE of predicted values at the gridpoints
%		in cross-validation.
% REss (nt x 1)r  reduction of error statistic from cross-validation of
%		the SSR models.  REss is a verification statistic.  
% RSQ1 (nmods x qq)r  R-squared for predictand PC scores from PCR
%		models.  Rows correspond to spatial models, columns to 
%		retained predictand PC's. RSQ1 is a calibration statistic.
% RSQ2 (nmods x q)r   "quasi" R-squared for spatial reconstruction
%		at each gridpoint or region.  RSQ2 is computed from the 
%		observed values and errors at the gridpoints by :
%
%			RSQ2 = 1.0 - (var(error) ./ var(observed));
%
%		RSQ2 could concievably be negative for some gridpoints.
%		RSQ2 is a calibration statistic.
% RSQ3 (nmods x q)r  squared correlation coefficient between
%		predicted and actual climate for each gridpoint or region.
%		RSQ3 is a verification statistic computed from the 
%		agreement of actual and observed climate for the 
%		cross-validation estimates.
% STATS (nt x 4)r  calibration statistics summarizing SSR models.
%		Columns: R-squared, adjusted R-squared, overall-F for
%		the regression equation; p-value of the overall-F
% YO (mTC x q)r  observed regional climate variable for the 
%		calibration period at the q gridpoints or regions.  YO is
%		either taken directly from input matrix YY, or computed
%		as described elsewhere by averaging the station climate
%		variable over stations.  The predictands for the PCR
%		models are PC scores of the data matrix YO.
% YP (mrec x q)r reconstructed regional climate variable at the
%		q gridpoints or regions for the mrec reconstruction years.
%		The vector yr4 holds the corresponding years.
% filesv (1 x ?)s  name of .mat file holding the output
% hf (mrec x 1)L  flag for extrapolations in PCR reconstructions;
%		extrapolation is flagged by "1", other values are "0"
% kcv (1 x 1)s  Y or N:  the output saved was generated with
%		cross-validation option on (Y) or off (N)
% kopt (1 x 4)i -- program options *** SEE INPUT
% lg (3 x 1)i -- number of delays, negative lags, positive
%		lags for SSR models  *** SEE INPUT
% mineigs (1 x 1)i minimum number of predictor PC's to be used
%		in PCR models *** SEE INPUT 
% modnum (mrec x 1)i  time series vector indicating which PCR model
%		applies for each reconstruction year;  go to the indicated
%		row of L6 to see which SSR's apply to a model;  go to the
%		indicated row of P2 to see year range for a particular model
% modyrs (nmods x 1)i  number of reconstruction years each PCR model
%		applies to
% nmods (1 x 1)i  number of PCR models
% tcrit (1 x 1)r  critical significance level (e.g., .95) for t-test
%		for entry of PC scores as predictors in PCR models 
%		*** SEE INPUT
% yr3 (mTC x 1)i years vector for calibration period
% yr4 (mrec x 1)i years vector for reconstruction period
%

%*********** CROSS-REFERENCE OF SAVED VARIABLES WITH OPERATIONS
%
%+++++++ SINGLE-SITE RECONSTRUCTIONS (SSR'S), CALIBRATION 
%
% ICD --index pointer: grouping of climate stations with trees
% DS -- observed tree-weighted local climate time series,calibration pd
% DHS -- single-site reconstructed time series, calibration period
% DHL -- (not currently saved in output .mat file) 
%		long-term single-site reconstructions, same size as input T
% Pyrs -- start and end years for single-site reconstructions;
% Fcrit --specified critical significance level for overall-F
% STATS -- R-squared, adjusted R-sqd, overall-F, p-value for F
% Cwgt -- regression coefficients for models
% IDX -- index pointer: averaging over SSR's
%
%+++++++ SINGLE-SITE RECONSTRUCTIONS (SSR'S), CROSS-VALIDATION 
%
% Mcv1 -- calibration means of the SSR's
% DHkey -- SSR single-site cross-validation estimates for 
% R2 -- calib R-squared for the mTC leave-5-out SSR models
% R2ss -- squared correlation coefficient, actual and observed
%		and predicted 
% REss --  reduction of error statistic   
%
%+++++++ SPATIAL MODELS (PCR MODELS), CALIBRATION 
%
% ICY -- grouping of climate stations to regions
% YO --  observed regional climate series, calibration period
%  yr3 -- corresponding year vector
% YP --  reconstructed time series, mrec reconstruction years
%  yr4 -- corresponding year vector
% modnum -- time series of spatial models
%   L6 -- associated variables in spatial models
%   nmods -- number of spatial models
%   mineigs -- minimum number of predictor PC's to be used
%   tcrit --  critical significance level for entry of predictors
%	 modyrs -- number of years in each spatial model (cv)
% EV -- "explained variance" statistic  
% NV -- summary information on predictors and predictands
% P2 -- start and end reconstructed year
% RSQ1 --  R-squared , based on PC-score predictands
% RSQ2 --  "quasi" R-squared, gridpoint
% hf -- flag for extrapolations
%
%+++++++ SPATIAL MODELS (PCR MODELS), CROSS-VALIDATION
%
% EVcv -- Like EV, but for calibration of each leave-5-out model
% PEbar -- root-mean-square error of prediction 
% RE -- reduction-of-error statistic
% RSQ3 (nmods x q)r  squared correlation coefficient between
%
%
% 
%********* INTERESTING VARIABLES NOT SAVED IN OUTPUT .MAT FILE
% To access these variables, search for "KEYBOARD-1" string, 
% de-comment the "keyboard" command, run the program, and grab
% the variables before going on.
%
%
% ns1 -- cv, number of clim stations averaged around each tree
% ns2 -- cv, number of clim stations averaged in each region
%		note:  ns2==[] of kopt(4)==2
% yr1 -- cv, year vector, input time series T, long single-site
%		reconstruction matrix DHL
% yr2 -- cv, year vector for input matrix C
% yr6 -- cv, year vector for optional input matrix YY
% DHL -- long-term matrix of single-site reconstructions
%
%
%
%*******************  GRAPHICAL SUMMARY OF RESULTS ************
%
% To use the following summaries, run spatrec1.m, save the results
% in a .mat file as prompted, clear the workspace, and load the
% .mat file.  You then have access to the required matrices.
%
%                        ****
% tsp plot of reconstruction with RMS pred error bars; wo/extrap flags
%
% e=pltrec1(yr3,YO(:,1),yr4,YP(:,1),PEbar(:,1),modnum,hf,ktrans)
%    where ktrans==1 for transparency or 0 for black bkgd screen plot
%
%                         ****
% tsp plot of reconstruction with RMS pred error bars;w/extrap;
%		screen prompted options for converting units (e.g, in to cm)
%		and truncating start year of plot
%
% e=pltrec2(yr3,YO(:,1),yr4,YP(:,1),PEbar(:,1),modnum,hf)
%
%                         ****
% tsp's summarizing change of predictor subsets and reconstruction
% 		accuracy with time.  Two sets of axes appear in one figure
%		window.  Top:  R2,RE (top) 
%			Bottom: % number of trees available, number of tree-ring 
%				PC's in pool of potential predictors, and 
%				number of PC's used in reconstruction model
%
% e=pltrec3(yr4,modnum,EV,RE,NV,ktrans)
%    where ktrans is defined as for pltrec1.m
%
% 
%                        ****
% Map the spatial distribution of errors. For this, would need 
% x,y coordinates for the predictand regional climate variables.
% A series of maps summarizing spatial variability is then possible:
%
% RSQ2 = 1 - (var(errors)/var(actual)); Calibration-period accuracy
%		with a statistic similar to R-squared.  Note, however, that
%		the regression R-squared applies to PCR models with PC scores
%		as predictands, while RSQ2 is in terms of the original variables
%		at their q gridpoints or regional points
% 
% RE reduction of error statistic at the q points.  This field
%		can be compared with RSQ2 for a summary of calibration and 
%		validation accuracy;  other cross-validation maps possible are
%
%		RSQ3 -- squared correlation coeff between actual and predicted
% 		PEbar -- square root of the average squared prediction error
%
%
%
%
% *****************  TABLE SUMMARY OF RESULTS **************
%
% The script files described here produce summary tables in tab-
% delimited form that can be easily imported to MS Word to make
% Footnote MS-Word files are available for slapping on the end
% of the MS-Word tables.  To use the script files, clear the 
% workspace, load the .mat file with the model results, and run
% the scripts.  You will be screen prompted for file info.
%
%                    ****
% spatout1.m  ...summary statistics on single-site reconstructions
%
%                    ****
% spatout2.m ... summary statistics of spatial reconstructions


%****************** STEP-1: ARGUMENT CHECK, SIZING, ALLOCATION


[mtemp,ntemp]=size(kopt);
if mtemp~=4 | ntemp~=1
	error('kopt must be 4 x 1 vector')
end
if ~(kopt(4) ==1 | kopt(4)==2),
	error('kopt(4) must be 1 or 2')
end


[n,nt]=size(T);
[mc,nc]=size(C);
[m2,n2]=size(YRS);
[m1,n1]=size(ICD);
[m4,n4]=size(IDX);
if kopt(4)==2, % specify regl climate matrix rather than compute from C
	if ~exist('YY'),
		error('YY does not exist, yet you set kopt(4)==2')
	end
	[m6,n6]=size(YY);
	if m2~=5| n2~=2, error('YRS must be 5 x 2'), end
	if ~exist('ICY'); % need to have a dummy matrix for storage
		ICY=[];
	end
else
	[m3,n3]=size(ICY);
	if m2~=4| n2~=2, error('YRS must be 4 x 2'), end
end

% delay and number of neg, pos lags for predictors are in lg
[mtemp,ntemp]=size(lg);
if mtemp~=3 | ntemp~=1, error('lg must be 3 x 1'), end


% each tree should be represented by a column of T and row of ICD
if nt ~= m1, error('col size of T must equal row size of ICD'), end


if YRS(3,2) > YRS(1,2),
	error('End year of calib pd later than end year of tree rings')
end
if YRS(1,1)>YRS(3,1)-lg(2),
	error('Start year of T too late for specified calib pd/ lags')
end
if YRS(1,2)<YRS(3,2)+lg(3),
	error('End year of T too early for specified calib pd/ lags')
end


if YRS(2,1)>YRS(3,1)
	error('Start year of C too late for specified calib pd')
end
if YRS(2,2)<YRS(3,2),
	error('End year of C too early for specified calib pd')
end
if  (YRS(4,1)-lg(2))<YRS(1,1),
	error('Recon period (lag-adjusted) too early first year')
end
if YRS(4,2)>= YRS(3,1),
	error('Recon period overlaps calibration period')
end
if kopt(4)==2; % YRS has row 5 with years of input regional clim mtx
	if YRS(5,1)>YRS(3,1)
		error('Start year of YY too late for calib pd')
	end
	if YRS(5,2)<YRS(3,2)
		error('End year of YY too early for calib pd')
	end
end


% Some key year-index pointers
yr1 = (YRS(1,1):YRS(1,2))'; % years vector for T
yr2 = (YRS(2,1):YRS(2,2))'; % years vector for C
yr3 = (YRS(3,1):YRS(3,2))'; % years vector for calib period
yr4 = (YRS(4,1):YRS(4,2))'; % years vector for recon period
if kopt(4)==2,
	yr6 = (YRS(5,1):YRS(5,2))'; 
end

L1T = yr1>=YRS(3,1) & yr1<=YRS(3,2); % to cal pd years in T
L1C = yr2>=YRS(3,1) & yr2<=YRS(3,2); % to cal pd years in C
L2T = yr1>=YRS(4,1) & yr1<=YRS(4,2); % to recon pd years in T
if kopt(4)==2,
	L1YY = yr6>=YRS(3,1) & yr6<=YRS(3,2)% to cal pd years in YY
end

%****************** STEP-2: MAKE SITE-AVE AND REGIONAL-AVE CLIM SERIES


%*** AVERAGE CLIMATE SERIES AROUND TREE SITES

disp('Starting to average sites into groups, stns into regions')
% Number of stations to average around each tree
L1=~isnan(ICD);
if n1==1; % special case, ICD is not a matrix but a cv
	ns1=L1;
else
	ns1=   (sum((L1)'))'; % col vector of number of clim station each tree
end
if any(ns1<1), error('Some row of ICD is all NaN'), end

% Form logical matrix corresponding to ICD that can be used in 
% matrix multiplication 
L2 = ind2log(ICD,nc);
D = C * L2';  % tsm of sum of ppt over stations near each site

% Check for NaNs in selected series
if any(any(isnan(D(L1C,:)))),
	error('NaN found in calib part of single-tree climate matrix D')
end

temp1=ns1'; % convert to row vector -- number of stations each site
denom=temp1(ones(mc,1),:); % expand to matrix
D = D  ./ denom;  % climate variable averaged over specified
		% stations near each tree;  D is the tree-centered climate
		% series to be used as the predictand variable in single-
		% tree reconstructions


% *** "GET" A REGIONAL CLIMATE MATRIX TO BE USED AS PREDICTANDS IN 
% RECONSTRUCTION.  THIS IS DONE IN ONE OF TWO POSSIBLE WAYS:
% (1) BY AVERAGING THE LOCAL CLIMATE VARIABLE OVER STATIONS, OR
% (2) BY READING IN THE REGIONAL CLIMATE VARIABLES.

% METHOD 1: AVERAGE THE LOCAL CLIMATE VARIABLE OVER STATIONS

if kopt(4)==1,  % This section applies conditionally
% Compute ns2, a column vector specifying how many stations to 
% average into each regional climate series
L1=~isnan(ICY);

% Handle special case of ICY a cv rather than a matrix.  This means 
% there is only one climate station per "region", and L1 is
% a column vector of ones
if n3 == 1; 
	ns2=L1;
else
	ns2=   (sum((L1)'))'; % col vector of number of clim station each tree
end
if any(ns2<1), error('Some row of ICY is all NaN'), end


% Form a logical matrix corresponding to ICY that can be used in 
% matrix multiplication 
L2 = ind2log(ICY,nc); % Note that nc is col size of C
Y = C * L2';  % tsm of sum of local climate variable over stations
%		in each region

% Check for NaNs in selected series
if any(any(isnan(Y))),
	error('NaN found in regional climate matrix Y')
end

% Divide the regional-sum precipitation by the number of stations
% in each region to get the regional average
temp1=ns2'; % convert to row vector -- number of stations each region
denom=temp1(ones(mc,1),:); % expand to matrix
Y = Y  ./ denom;  % climate variable averaged over specified
		% stations in region;  Y holds the regional-average climate
		% series to be used (indirectly) as the predictand variables in 
		% spatial recontruction. Note that Y(L1C,:) will be the calib-pd
		% rows of Y

end ; % of METHOD 1 for regional climate matrix


% METHOD 2: USE INPUT MATRIX YY AS REGIONAL CLIMATE SERIES

if kopt(4)==2,
ns2=[];
Y=YY;   % Note that Y(L1YY,:) holds the calib-period rows
		% of Y
end; % of "METHOD 2" if for regional climate series formation


[my,q]=size(Y);  % spatial model will have q predictands
disp('Finished averaging into groups and regions')

%********************  STEP-3:SINGLE-SITE RECONSTRUCTIONS  *************
%
% Check validity of time periods and lags
% Form lagged-predictor matrix
% Initialize matrices
% Estimate parameters of models
% Reconstruct


disp('Starting single-site reconstructions')

% Check validity of time periods and lags 
% Spatrec1.m is only for t-1,t,t+1 models
if lg(1)>0, error('delay lg(1) must be zero or negative'), end
if lg(2)<0 | lg(3)<0,
	error ('lg(2) and lg(3) must be non-negative')
end
if lg(2)~=1 | lg(3)~=1,
	error ('lg(2) and lg(3) must equal 1 for this program')
end


% Form lagged-predictor matrix;  TL will have lag-0 data in
% first cols, followed by data for lags -1, then for +1
% Spatrec1 handles only t-1,t,t+1 models

[TL,yrTL] = lagyr3(T,YRS(1,:),lg);
TC = TL(L1T,:);  % subset of TL to be used for full calibration
[mTC,nTC]=size(TC);

% Initialize matrices
% Allocated space for several of these matrices is set by 
% t-1,t,t+1 restriction of spatrec1.m
a=NaN;
Cwgt=zeros(4,nt); % will hold regresion coefs
STATS = a(ones(nt,1),ones(4,1)); % calibration stats of SSR models
DS = a(ones(mTC,1),ones(nt,1)); % calib-pd observed predictands
DHS = a(ones(mTC,1),ones(nt,1)); % calib-pd recons
DHL = a(ones(n,1),ones(nt,1)); % entire-rows-of-T recons
II1 = a(ones(nt,1),ones(3,1)); % Cols of TL for each SSR
II2 = a(ones(nt,1),ones(3,1)); % Relative columns of selected
		% predictors for SSR models

% Estimate parameters
% Loop over each tree
for i = 1:nt;
	c1 = zeros(4,1); % To hold regr coefficients for this tree
	I1 = [i  i+nt   i+2*nt]; % for lags 0, -1, +1
	II1(i,:) = I1;
	y = D(L1C,i); % predictand
	DS(:,i)=y;  % store calib-period actual local-climate series
	if any(any(isnan(TC(:,I1))))
		error('A calib-period predictor value is NaN')
	end
	[I2,I4,stats,c,e,yhs]=stepr2(TC,y,I1);
	
	II2(i,1:length(I2))=I2;
	lenc = length(c);
	c1(1)=c(1);
	c1(I4+1)=c(2:lenc);
	Cwgt(:,i)=c1;

	% Is P-value for computed F sufficiently low?
	Fcompute=stats(3);
	Fp=1-fcdf(Fcompute,(length(c)-1),(mTC-length(c)));
	if Fp>1-Fcrit,
		disp(['Poor overall-F; pvalue > threshold; tree ',int2str(i)])
		pause(3)
	end
	
	STATS(i,:) = [stats Fp];
	DHS(:,i) = yhs;
	Tall = [ones(n,1) TL(:,I1)];
	DHL(:,i) = Tall * c1; 
end

% DHS now holds the single-site recons for the calib period
% DHL now holds recons for all years covered by T




Pyrs1=yrsgosp(DHL(L2T,:),YRS(1,1)); % start,end years of non-NaN data for
% SSR's, considering pre-calibration period only
Pyrs2=yrsgosp(DHL(L1T,:),YRS(3,1)); % Likewise, but considering only the
% calibration period

Pyrs=[Pyrs1(:,1) Pyrs2(:,2)];  % start,end years on any non-NaN data in DHL
disp('Finished single-site reconstructions')


%****** STEP-4   SPATIALLY AVERAGE SSR'S OVER TREES

disp('Beginning To Spatially Average Single-Site Reconstructions')

L3=~isnan(IDX);
if n4 == 1; % special case: IDX Is not a matrix but a cv
	ns3=L3;
else
	ns3=   (sum((L3)'))'; % col vector of number of sites to average
end
if any(ns3<1), error('Some row of IDX is all NaN'), end

% Form logical matrix corresponding to IDX that can be used in 
% matrix multiplication 
L4 = ind2log(IDX,nt);
LX8 = isnan(DHS);
LX9 = isnan(DHL);
sum8 = sum(sum(LX8));
sum9=sum(sum(LX9));
DHSz=DHS;
DHLz = DHL;

% If any NaNs in DHSz or DHLz, replace the NaNs with zeros
% This will let averaging over trees be possible using
% matrix multiplication possible
if sum8>0;
	DHSz(LX8)=zeros(sum8,1);
end
if sum9>0;
	DHLz(LX9) = zeros(sum9,1);
end


X = DHSz * L4';  % tsm of sum of single-tree recons over sites
Xbig = DHLz * L4';
% Xbig has values or zeros as elements.  Some of the values might
% not be averages over all the specified sites as indicated by
% IDX. Must check this and replace invalid elements of Xbig with
% NaN

% Make a matrix like DHL, but with non-NaN elements filled with
% ones and all other elements zeros
DHLvalid = ~isnan (DHL);

% Make a matrix same row-size as DHL and DHLvalid, with 
% sum of number of sites in region as elements.  All rows same.
L4sum=sum(L4');
L4sum=L4sum(ones(n,1),:);

B1 = DHLvalid * L4';
LL4 = B1 ~= L4sum;

sum12 = sum(sum(LL4));

Xbig(LL4) = a(ones(sum12,1),:);



% This obsolete
%a = NaN;
%if sum8~=0;
%	 X(LX8) = a(ones(sum8,1),:);
%end
%if sum9~=0;
%	Xbig(LX9) = a(ones(sum9,1),:);
%end

% Check for NaNs in X
if any(any(isnan(X))),
	error('NaN found in predictor matrix X')
end


% Check that Xbig is not all NaN for the rows marked by YRS(4,1)
% as start of reconstruction period
L15 = (~all(isnan(Xbig')));
i15 = find(L15);
yrgo = yr1(i15(1));
if yrgo>yr4(1);
	disp('First year with any valid data as averaged over series')
	disp([' as marked in IDX is ',int2str(yrgo)]);
	disp([' But you have specified ',int2str(YRS(4,1)),' as the '])
	disp(' first year to be reconstructed');
	error(['Revise YRS(4,1) to ',int2str(yrgo)]);
end

temp1=ns3'; % convert to row vector -- number of sites in each group
denom = temp1(ones(mTC,1),:); % expand to matrix
X = X ./ denom;

denom=temp1(ones(n,1),:); % expand to matrix
Xbig = Xbig  ./ denom;
disp('Finished spatially averaging single-site recons')

%********** STEP-5: IDENTIFY SSR SUBSETS *****************
%
%******* Get col-index matrices for distinct recon-period data sets

disp('Starting to identify SSR subsets')
YRS1 = [YRS(1,:); YRS(4,:); YRS(3,:)];
[L6,modnum]=uniqmods(Xbig,YRS1);
% L6 (mL6 x nL6)L pointer matrix telling which cols of X and Xbig
%		are active for each of mL6 models
% modnum (matrix)i which model (row of L6) applies for each
%		of the reconstruction years in Z as pointed to by 
%		row 2 of yrs
[nmods,dum1]=size(L6);  % nmods is number of different models
disp('Finished identifying SSR subsets')


%******** STEP-6:  SPATIAL (PCR) RECONSTRUCTIONS ********

disp('Starting PCR modeling')
% Loop over each of the recon-year models, estimating parameters,
% getting reconstructed values, and saving model statistics

% Allocate 
anan=NaN;
mrec = sum(L2T); % number of reconstruction-period years
YP = a(ones(mrec,1),ones(q,1)); % reconstructed climate
hf = a(ones(mrec,1),:); % extrapolation flag
EV=anan(ones(nmods,1),:);
RSQ1=anan(ones(nmods,1),ones(q,1));
RSQ2=anan(ones(nmods,1),ones(q,1));
NV=anan(ones(nmods,1),ones(5,1));
modyrs=anan(ones(nmods,1),:); % number of years in each
%		spatial model

% Get the calibration-period predictand data, "Y Observed"
if kopt(4)==1,
	YO = Y(L1C,:);
else
	YO = YY(L1YY,:);
end

for n5=1:nmods;% 2;
	disp(['Beginning spatial model ',int2str(n5)])
	jcols =L6(n5,:);
	% Find out rows (years) in T to be reconstructed with this model
	L7 = modnum==n5;
	modyrs(n5)=sum(L7);
	XXbig = Xbig(L2T,:);
	Xrec = XXbig(L7,jcols);
	Xcal = X(:,jcols);

	[yp,yrec,rsq1,rsq2,ev,nv,h1]=osr(Xcal,YO,Xrec,mineigs,tcrit,kopt);
	hf(L7)=h1; % extrapolation flag
	NV(n5,:)=nv;
	RSQ1(n5,1:nv(2))= rsq1;
	RSQ2(n5,:)=rsq2;
	EV(n5) = ev;
	YP(L7,:) = yrec; 
end

% clear up some memory
clear T Xbig

% Make 2-column matrix with first and last year for each spatial 
% model.
AA = modnum(:,ones(nmods,1));  % expand cv modnum to matrix by
		% duping rows
BB = (1:nmods); % rv 
BB = BB(ones(length(modnum),1),:); % expand rv to matrix
BB=AA~=BB;
nAB = sum(sum(BB));
aa=anan(ones(nAB,1),:);


AA(BB)=aa; % AA now has NaNs except where 
		% spatial models apply

P2 = yrsgosp(AA,yr4(1));

clear AA BB nAB

disp('Finished the spatial (PCR) reconstructions -- and calibration')

%keyboard;  %KEYBOARD-1

%***********************  CALIBRATION COMPLETED ******************

kcv=[];
quest1=questdlg('Continue with crossvalidation?');
switch quest1;
	case 'No';
		kcv='N';
	case 'Cancel';
		kcv='N';
	case 'Yes';
		kcv='Y';
end


if kcv=='Y';

disp('Beginning Crossvalidation')



%********  STEP 7: CROSSVALIDATION -- SINGLE TREE  **********
%
% Must repeat the estimation of single-tree models using the same
% models as in the previous single-site reconstruction, but
% leaving out 5 years from the estimation. the 5 years are
% successively shifted by one and the model recalibrated
%
% Have nt different trees for which single-tree regressions are done.
% Regression models were previously estimated for these trees. We will
% use the same designated predictor variables as found in the full
% calibration.  For example, if the full calibration indicated a
% no-lag single-site reconstruction model for tree#1, a no-lag model
% will be assumed for the leave-5-out models.  Of course, the 
% regression parameters will be re-estimated.
%
% Outer loop over the the mTC leave-5-out calibration periods.
% Inner over the nmods spatial models.
% 
% Within the outer loop, must re-calibrate the single-tree-models.
%
% Critical  output info includes, for each of the nummods models
% - mean, median, lo, hi R-squared for the leave-5-out calibrations
%		in single tree mods
% - r-squared and RE for the single-tree validations
% - mean EV for the leave-5-out spatial regressions
% - validation estimates of climate at the q points
% - RV computed from validation-period residuals
% 
%
% Keep these settings from the earlier steps
% - same predictor variables for the single-tree regressions
% - same number of eigenvectors retained on the predictands and
%    predictors (qq and pp from pcareg1.m).  Note that the
%    exact predictors of the pp are not set in stone. This seems
%    reasonable, because the same eigenvectors might not result
%    from PCA on the full calib period and on the full period
%    minus the years needed for the validation estimate
%
%  Currently available variables from earlier steps:
%
% TC (mTC x nTC)r lagged tree-ring data, mTC calib years,
%		nt*3=nTC cols (lags 0,-1,+1)
% II2 (nt x 3)i  cols of TC for models, NaN right-filled
% NV (nmods x 5)i [q qq p pp nused] number of predictands, retained
%	predictands, predictors, retained predictors, and
%	predictors used x nummods)i 

% Initialize key matrices
R2 = a(ones(nt,1),ones(mTC,1)); % Calib R-squared for the mTC
	% leave-5-out single-site models at each of the nt trees
EVcv=a(ones(mTC,1),ones(nmods,1)); % calibration stat
DHkey=a(ones(mTC,1),ones(nt,1));% single-site predictions from 
	% crossvalidation, for the period outside calib period
Z = a(ones(nmods*mTC,1),ones(q,1)); % predicted climate at q points 
zgo = 1:mTC:(nmods*mTC-(mTC-1)); % starting row in Z of each 
	% mTC-row predicted submatrix
Mcv1=a(ones(mTC,1),ones(nt,1)); % Calibration-period means of single-
	% tree local climate data used for models that produce the 
	% crossvalidation estimates.  This matrix needed for computing
	% reduction of error statistic (see RE1) for single-tree models
Mcv=Z; % calibration-period means needed for RE statistic


% Get the index for crossvalidation years
LTC=crospul2(mTC,1,1);


for k1 = 1:mTC;  % Loop over crossvalidation periods
	disp(['Starting crossv period ',int2str(k1)]);
	Lr = LTC(:,k1); % rows of calibration period to calibrate model on 
	yr5 = yr2(Lr); % year-cv for current leave-5-out calib period
	nobs = sum(Lr); % number of calibration years
	TC1 = TC(Lr,:); % predictor calibration data (lagged tree matrix)

	DHScv=a(ones(nobs,1),ones(nt,1)); % single-site predictions, 
		% cross-validation periods, calibration data

	% Single-site regression, looping over sites
	for k2 = 1:nt;
		y = D(L1C,k2); % full calib-period predictand vector
		y = y(Lr);  % just the sub-period for crossvalid. fitting
		Mcv1(k1,k2)=mean(y); % predictand mean for this calib pd and
			% site 
		ii2 = II2(k2,:);
		npreds = sum(~isnan(ii2));
		ii2 = ii2(1:npreds);
		[bwt,r2,yhat] = regress1(y,[ones(nobs,1) TC1(:,ii2)]);
		 
		Tkey = [1   TC(k1,ii2)];
		DHkey(k1,k2) = Tkey * bwt; % prediction for the indep data
		R2(k2,k1)=r2; 
		%Tip: plot(R2); shows variability of R-squared from one
		% model to another as different groups of obs are omitted
		% [STATS(:,1) (R2(:,k))] shows comparison of calib R-squared for
		%		full-period calibation with calib-R-squared for
		%		kth leave-5-out calibation
		DHScv(:,k2) = yhat; % calib-period predictions for currrent
			% leave-5-out model
		% Tip: For time series plots of (1) observed single-site
		% climate series, (2) full-period calibration-period 
		% predictions, and (3) leave-5-out calibration-period
		% predictions for current model, do this:
		%  plot(yr2,D(:,3),yr2,DHS(:,3),yr5,DHScv(:,3))
	 end

	% Get crossvalidation statistics for single-site models
	% R2ss is the squared correlation coefficient between observed
	% and predicted local climate variable.  REss is the reduction of
	% error statistic.  R2ss and REss are column vectors of length 
	% equal number of sites.
	[R2ss,REss]=cvstat1(DS,DHkey,Mcv1);

	% Spatially average single-tree reconstructions
	% Form logical matrix corresponding to IDX that can be used in 
	% matrix multiplication -- L4 formed previously
	LX10 = isnan(DHScv);
	dhkey = DHkey(k1,:);
	LX11 = isnan(dhkey);
	sum10 = sum(sum(LX10));
	sum11 = sum(LX11);
	DHScvz=DHScv;
	dhkeyz=dhkey;

	
	if sum10~=0,
		 DHScvz(LX10)=zeros(sum10,1);
	end
	if sum11~=0;
		 dhkeyz(LX11)=zeros(sum11,1);
	end

	% Weight single-site predictions over sites
	Xcv = DHScvz * L4';  % tsm of sum of single-tree recons over sites
	Xkey = dhkeyz * L4';

	a = NaN;
	if sum10~=0
		Xcv(LX10) = a(ones(sum10,1),:);
	end
	if sum11~=0,
		Xkey(LX11) = a(ones(sum11,1),:);
	end

	temp2=ns3'; % convert to row vector -- number of sites in each group
	denom = temp2(ones(nobs,1),:); % expand to matrix
	Xcv = Xcv ./ denom;
	Xkey = Xkey ./ temp2;

	% Check for NaNs in selected series
	if any(any(isnan(Xkey))),
		error('NaN found in predictor matrix Xkey')
	end

	if any(any(isnan(Xcv))),
		error('NaN found in predictor matrix Xcv')
	end

	% PCA REGRESSION FOR CROSSVALIDATION	
	 
	for n5=1:nmods; % Loop over the different tree-group models
		disp(['   Starting tree-group model ',int2str(n5)]);
		jcols =L6(n5,:);
		Tcal = Xcv(:,jcols);
		Ycal = YO(Lr,:);
		Trec = Xkey(jcols);
		[ypcv,z1,rsq1,rsq2,ev,nv]=osr(Tcal,Ycal,Trec,mineigs,tcrit,kopt);
		%RSQ1cv(n5,1:nv(2))= rsq1;
		%RSQ2cv(n5,:)=rsq2;
		EVcv(k1,n5) = ev;

		 iz = zgo(n5) + k1-1; % row index in Z to store this prediction
		Z(iz,:) = z1; 
		Mcv(iz,:) = mean(Ycal); % calibration-period means needed for RE
  end

end;  % k1 loop over crossvalidation periods

% $$$ "keyboard" special $$$
% To compare time series of cross-validation estimates for 
% single-site reconstructions with the observed series, do the
% following in a loop over k
%[(DHkey(k,:)  D(k,:)]  ... to view data
% plot(yr2,DHkey(k,:),yr2,D(k,:)) ... to check agreement with plot


% Compute verification: squared correlation coefficient and RE
% Initialize matrices
RE = a(ones(nmods,1),ones(q,1)); % reduction of error at each gridpoint
RSQ3=zeros(size(RE));  % squared correl coeff between predicted and actual
	% climate for each spatial-model/gridpont
MSE=a(:,ones(q,1)); % initialize
MSM=zeros(size(MSE));
PEbar=zeros(size(RE)); % sensible estimate of prediction error

% Get calibration-period observed regional climate series
if kopt(4)==1,
	YO = Y(L1C,:);
else
	YO = YY(L1YY,:);
end


for i= 1:nmods ; % loop over spatial models
	Y1 = Z(zgo(i):(zgo(i)+mTC-1),:); % predicted climate
	M1 = Mcv(zgo(i):(zgo(i)+mTC-1),:); % calib-period means
	EM = YO - M1; % actual minus calib-period mean
	EY = YO - Y1; % residuals = observed minus predicted
	PRESS=sum(EY .* EY); % prediction error sum of squares, rv
	PEbar(i,:) = sqrt(PRESS/mTC); % Weisberg, top, p. 230 -- a "sensible
		% estimate of average prediction error"
	MSE = mean (EY .* EY); % rv of mean square errors at the 
	MSM = mean(EM .* EM);

	% RE statistic
	RE(i,:) = 1 - (MSE ./ MSM);

 	% Correlation coefficients
	for j = 1:q;
		yreal=YO(:,j);
		y1 = Y1(:,j);
		rr = corrcoef([y1 yreal]);
		RSQ3(i,j) = rr(1,2) .^2; 
	end
	
end

end ; % of if kcv



%******************  SAVE OUTPUT IN A MAT FILE  ***************

quest1=questdlg('Save output in a .mat file?');
switch quest1;
	case 'No';
		ksv='N';
	case 'Cancel';
		ksv='N';
	case 'Yes';
		ksv='Y';
end

if ksv=='Y' ;  % You want to save key output 
	SET1=[' save scratch1 DS DHS STATS Fcrit Cwgt Pyrs'];
		% from single-tree regressions
	SET2=[' save scratch2 EV L6 NV RSQ1 RSQ2 YO YP kcv '];
	SET3=[' save scratch3  modnum nmods modyrs P2'];
	SET4=[' save scratch4  ICY IDX ICD mineigs tcrit lg kopt '];
	SET5=['save scratch5  yr3 yr4 '];
	SET6=[' save scratch6 RE RSQ3 hf PEbar EVcv DHkey Mcv1 R2ss REss R2'];

	

	eval(SET1)
	eval(SET2)
	eval(SET3)
	eval(SET4)
	eval(SET5)
	if kcv=='Y';
		eval(SET6)
	end

	clear
	
	[file4,path4]=uiputfile('*.mat','Output .mat file to store results in');
	pf4=[path4 file4];


	load scratch1
	load scratch2
	load scratch3
	load scratch4
	load scratch5
	if kcv=='Y';
		load scratch6
	end
	
	eval(['save ' pf4]);
else
end 

% Clean up by deleting scratch files
dos('del scratch1.mat');
dos('del scratch2.mat');
dos('del scratch3.mat');
dos('del scratch4.mat');
dos('del scratch5.mat');
if kcv =='Y',
	dos('del scratch6.mat')
end




% prew1.m   estimates parameters for specified arma model.  
%	Whitens time series using the estimated arma model.
%	Designed to use input from arma.m, which does exploratory
%	model fitting for low-order arma models.

% 8-7-90, by D. Meko, Laboratory of Tree-Ring Research, U of AZ 85721

%###########  DESCRIPTION   #######################################
%
% Whitens an array of time series, using low-order arma models that
% have previously been selected using arma.m.  As with arma.m, this
% program was designed to consider as candidate models only  (1,0)
% (2,0), (0,1), and (1,1).
%
% The output array has columns of original data replaced by the 
% whitened data.  If a "null" whitening model was selected, the
% respective cols in the whitened array are just the original data.
%
% See description of arma.m for references.
%
% The main steps are:
%
% 1. Read in array of time series, and arrays specifying models.
% 2. Estimate model parameters.
% 3. Fill appropriate cols of an output array with the arma residuals.
% 2. 

%####### PRELIMINARY SETUP ######################################
%
% 1. Run arma.m to select best low-order arma model.
% 2 Three arrays must be loaded before running prew1.m:
%	1. A - data array, with year as col 1
%	2. ends - ?x3 array , beginning and ending year for modeling
%	   		period for each time series.  Both A and ends are 
% 			also needed  (along with colwhich) to run arma.m.
%	3. mods - a ?x3 array specifying which arma model was selected
%			for each valid time series in A.  Array mods is created 
%			by arma.m.  If you leave MATLAB between running arma.m
%			and prew1.m, be sure to save mods.
% 2. Start MATLAB.
% 3. >load A, ends, and mods, if not already loaded.
% 4 >prew1... to start this program. 
% 
%
%###############  CLOSING ACTIONS  #############################
%
% 1. You may want to save one or more arrays;
%
%		E...residual, the whitened array
%		TBL1,TBL2... tables summarizing the fitted models
%		mods... the array instructing which models to fit
%		modstat1... an array with info on the selected models.
%			Note that mods and mdstat1 are ouput of arma.m, not
%			prew1.m.
%
%
%###############  FEATURES   ####################################
%
% * Screen Plots of resids superposed on original data.
% * Flexibility to fit zero-order (no prewhitening) to selected series
% * Flexibility to use partial time segment of input array for model
%   estimation.  This segment can vary from series to series.  (As
%   yet, whitening is done only for that segment used to fit the model.)
% * Output tables that summarize model, variance explained, Porte Manteau,
%   Parameter estimates and standard errors.
% * Adaptable for hard-copy plots by insertion of "keyboard" commands
%   into program text.
% * Convenient for large data sets:  originally used to whiten tsa,
%   1700-1979, holding tree-ring indices for 249 chronologies.

%#############   LIST OF VARIABLES  #####################################
% 
% A...matrix of original time series data.  Each row a year, each col
% 	beyond the first a variable. Col 1 is the year.
% e...	residuals of an estimated AR model 
% ends..(ncols x 3)  tells valid period for analysis for each series.
%   First col is series number.  Second col is beginning valid row of A.
%   Third col is ending valid row.
% evar... variance of residuals of model fit.  Computed over only
%    those observations beyond the first nparam
% E...  matrix of residuals for all series.  Any cols (variables)
%	not whitened will still hold original data from A.  Any portion of
%	a whitened variable outside the whitening period will also have
%	original data values in the repective rows for those years
% k... loop index for which series to analyze
% m...  used in various places to hold row size of an array
% mods...? x3, order model to be fit to each series.  Col 1 is
%	series #, col 2 is ar order,  col 3 is ma order.  If both
% 	ar and ma orders are zeros, means don't fit a model.
%	Then the original data of A would be retained in that col of E.
%	Mods is input for prew1.m.  Mods is originally created by running
%	arma.m, which allows selection of the models.
% mseries... number of rows in mods, from call to size
% mz... mean value of original time series z
% n...  used in various places to hold col size of an array
% ncols... number of columns in A
% nparam... totl # of arma params in a model
% nrows... number of rows in A
% nyears ... number of valid years to analyze in series
% opt1... control over period to be analyzed
%	Y- same period (all rows) for all  series
%	N- ends holds beginning and ending years
% standards... 1x3, standard errors of ar-1 ar-2 ma-1 params
% TBL1...matrix model results: series #, AR order, MA order, 
%	orig var, resid var, ratio of resid variance to orig variance
% TBL2...series #, series mean, first ar param, second ar param,
%	first ma param, se of first ar param, se of second ar parm,
%	se of first ma param.  Zeros if not applicable.  For the null model
%	(arma(0,0)), only the series number is non-zero.
% th....packed matrix of parameter estimation results:
%  row 1 var(e), sampling interval T, nu,na,nb,nc,nd,nf,nk
%  row 2 FPE, year, mo, date, hr, min, and command model gen by
%  row 3 vector of est params: A,B,C,D, and F
%  row 4 to 3 +n:  estimated param covariance mtx, where n is sum of
%		all orders (number of est params)
% varz... variance of z, computed only after deleting number of
%   initial years equal to selected AR order.  This deletion to give
%   better comparison with variance of residuals.
% yr...	vector of years for the time series, read from col 1 of A 	
% yr1..   vector of valid years for a particular time series
% z...	cv, a time series to be analyzed; raw and mean removed

%#################  END OF LIST ####################

% A holds in the series to be whitened.  First col is the 
% year, remainder of cols are the variables.  Last few cols
% may be dummy variables needed to fill out array to a size best
% input to arma, prew, and prew2.

yr = A(:,1);
[nrows,ncols] = size(A);   % how many rows and cols in A

home
clc

%################# Pre-allocate some arrays #####################

E=A;  % this array will eventually hold arma residuals
		% For now, fill with the original data of A
E(:,1)=yr;  % The residuals array should also have year as first col.


% ############### A different set of years may be analyzed for each series.

clc, home
opt1 = input('Same rows to be analyzed for all series? Y/N [Y]: ','s');
if isempty(opt1) 
	opt1='Y';
else
end
disp(opt1);

[mseries,n]=size(mods);  % mseries is the number of time series to be
		% treated.  This number may include those evaluated as order
		% 0,0, and hence with no model to be estimated.


TBL1=zeros(mseries,6);
TBL2=zeros(mseries,8);


% ###### Loop through for each series  ######

for k = 1:mseries;

nn = mods(k,2:3);

%  Optionally use different analysis period for each column.

if opt1=='Y'
	rowgo=1; rowstop=nrows; nyears=nrows;
else
	rowgo=ends(k,2)- yr(1) +1;
	rowstop=ends(k,3) - yr(1) +1;
	nyears=rowstop-rowgo+1;
end


% Fill proper year and data into z and yr1

z=A(rowgo:rowstop,k+1);
yr1=yr(rowgo:rowstop);
	

% subtract the mean

mz=mean(z);
z = dtrend(z);



%pause
nn=mods(k,2:3);
nparam = sum(nn)  % total number of ARMA parameters in model



% Compute variance of original time series

varz=(std(z(nparam+1:nyears)))^2;



if sum(nn)~=0  % only fit a model if ar and ma orders are
						% not both zero.


% estimate the parameters of the selected model
% Must call different mfunction depending on whether 0 ar or ma params

if nn(2)==0
	th=ar(z,nn(1));
else
	th = armax(z,nn);
end

% compute the acvf of the residuals of model, and plot acf

[e,r] = resid(z,th); 
text(0.5,0.9,['Series #',num2str(k)])
text(0.5,0.7,['Model: ',num2str(nn(1)), ' - ',num2str(nn(2))])
pause(2)
% Compute variance of residuals over obs minus the first nparams,

evar=(std(e(nparam+1:nyears)))^2;   % variance of residuals

% compute variance of whitened series, pct var
% of original retained in the whitened series, and fill the table
% describing model results.  Note that r(1) already holds variance of 
% residuals.

% replace a col of the zeros matrix E with residuals from model,
% after adding back in the series mean which had previously been 
% removed

E(rowgo:rowstop,k+1)=e+mz;

plot(yr1,z,yr1,e)  
title(['Time series # ',num2str(k),', original-red, whitened-green'])  
pause(2)


TBL1(k,:)=[k  nn  varz  r(1)   (evar/varz) ];

standards=[0 0 0];
if (nn(1)==1 & nn(2)==0) 
	standards=[sqrt(th(4,1)) 0 0];
elseif (nn(1)==1 & nn(2)==1)
	standards=[sqrt(th(4,1)) sqrt(th(5,2)) 0];
elseif (nn(1)==2 & nn(2)==0)
	standards=[sqrt(th(4,1)) sqrt(th(5,2)) 0];
elseif (nn(1)==0 & nn(2)==1)
	standards=[sqrt(th(4,1)) 0 0];
else
	% no action
end
	
TBL2(k,:)=[k  mz  th(3,1:3) standards] ;

else  % Here if both specified ar and ma order were zero
	
TBL1(k,:)=[k 0 0 varz varz 1.0];
TBL2(k,:)= [k mz 0  0  0  0  0  0];

end   % End of if stmt referring to null model
end % End loop through each series 

clc
home
disp('    Series     ARO         MAO   Var(z)    Var(e)       Pct     ' )
fprintf('\n')
disp(TBL1), pause
clc 
home
disp(' Series mean, followed by AR parameters, then their SEs')
fprintf('\n')
disp(TBL2), pause

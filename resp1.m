% resp1.m  

% Written by Dave Meko, January 1992.

% Application:  Screening of a tree-ring index against monthly
% precipitation and temperature at a station to find out 
% which months' climate signals are most important to
% tree rings.

% Parts of pgm are modified from resp.m

% Objectives; tri=tree-ring index, p=precip, t=temperature
%
% 1. r between  tri and individual monthly p,t
% 2. r between ar-whitened versions of tri and ar-whitened
%	versions of p,t
% 3. r between tree-ring index and each of 24 individual whitened
%	p,t series, after each of these has had linear correl with
%	other 23 clim vbls removed.

%*************  METHOD   *********************************

% Tree-ring index whitened by fitting variable-order model and 
% selecting best model by AIC.  Monthly t and p similarly whitened,
% but with models fit to them.  Models for tri, p, t all fit for 
% approx. same time period.  If AIC for 0-order model is lowest,
% original time series is accepted as whitened version.

% Removal of linear correl for part (3) above is accomplished by
% stepwise mult linear regression, with up to 23 variables allowed
% to enter until adjusted R-squared is maximum.


%*****************   STEPS   ******************

% 1. Get appropriate subsets of rows and cols of tree and climate data.
% 2. Compute whitened version of the tree-ring index.
% 3. Compute whitened version of the precip-temperature array
% 4. Compute uncorrelated-residual climate array
% 5  Compute r between the unwhitened tree-ring and unwhitened climate
% 6  Compute r for corresponding whitened series
% 7. Compute r between whitened tree-ring and de-correlated, whitened
%	climate array.
% 8.  Preview the various bar plots.
% 9.  Optionally produce ascii plots for grapher use.


%*************   PRELIMINARY STEPS  ****************************************

% Run FIL999.F  and REAR1.M  to get the ppt and temperature data into MATLAB
%  array form.  Each array must have 13 numeric-only cols.  These arrays
%  must not have header lines.  First col is the year;  remaining cols are
%  the monthly values for Jan through Dec.

% Program doesn't handle missing values.  MATLAB also won't 
% recognize arrays with blank elements.  So put some dummy value in for missing
% values in arrays beforehand.  The dummy value must be numeric--like
% 09999.

% Program will prompt for period on which you want to run the correlation
% analysis, which can be subset of years of entire arrays A and B.  Make
% sure all the data in the period for analysis is valid (i.e., does not
% contain special missing-value codes, like "99999").

% The tree-ring array is assumed to be at least 2 columns.  First column
% holds the year, succeeding cols the tree-ring indices for one or more
% cores or sites.  As for the PPT and temperature arrays, no alphnumerics
% are allowed, no header lines, and no blank cells.

%*************** PRELOADS  ****************************************

% X.MAT	monthly precip array and related files in matlab form
%	A	(m1 x 13) the precip array
%	B	(m1 x 13) the temp array
%	C   (m2 x n2) the tree-ring array
% 	endmo - ending month of growth year; e.g., 8=august
%  nhi - max ar order to try to fit

%***********   RESTRICTIONS   ************************************

% A, B must be of same size and cover same years
% 
% C  	can cover any period overlapping A and B, but must have at least
%	nhi more years at the start to allow for loss of data when 
%	whitening.  
%
% C 	cannot have missing data in the period specified for analysis
% 	

%***************   GRAPHICAL OUTPUT    *****************************

% 1.  (2 x 2) bar: TL: AR-orders for each of 12 precip series
%					BL: corresp pct variance in explained by model
%					TR: AR orders for temperature series
%					BR: corresp pct var expld
%
% 2.  (2 x 2) bar
%
% 		TL: r between tree and p, with 95% sign dashed line
%		BL: r between tree and t, etc
%		TR: r between whitened versions of series in TL
%		BR: r between whitened versions of series in BL
%
% 3.  (2 x 2) bar
%
%		TL: R-squared for decoupling regression model, ppt series
%		BL: Likewise for temperature
%		TR: Corresp F-ratio for values in TL
%		BR: Corresp F-ratio for Values in BL
%
%		Two stars or one star above bar: 99%, 95% significance
%
% 4.  (2 x 2) Like figure (2), but for decoupled versions of 
%		whitened series.


%****************   SELECTED VARIABLES  (EXCEPT OUTPUT ARRAYS)   *********************
%
% d	(scalar) which tree-ring series to analyze--in col d+1 of C
% k - miscellaneous screen input: 
%   - save ascii files for surfer?
% k1 - if loop on scaling A,B, creating D, E
% k2 - for loop for treating a given tree-ring series
% k3 - code for which type of array to generate.  This code number
%    is suffixed to ascii file name of file F?.dat.  So, if array type
%    is 
% L1 (lcv) marks approp rows of D,E for this tree-ring series analysis
% L2 (lcv) marks approp rows of C
% m1,n1  of ppt array A and temperature array B
% m2,n2  of tree-ring array C
% n3     of series designator


%********  LOADS, PRELOADS, AND PREALLOCATES   *****************


% Check that necessary input exists

t = [exist('A') exist('B')  exist('C')   ...
    exist('endmo')  exist('nhi')];

if t~= [1 1 1  1 1 ]
	disp('YOU FORGOT TO LOAD APPROPRIATE ?.MAT')
	keyboard;
end

[m1,n1]=size(A);
[mtemp,ntemp]=size(B);
[m2,n2]=size(C);
if (mtemp~=m1 | ntemp ~= n1)
	error('A AND B ARE NOT OF SAME SIZE!');
end

n3=1;             % Number of series to analyze, change later to while lp

for k2=1:n3;  % loop for each tree ring series to be analyzed

	d=input('WHICH TREE-RING SERIES? ');
	a(1)=input('FIRST YEAR FOR ANALYSIS PERIOD: ');
	a(2)=input('LAST YEAR FOR ANALYSIS PERIOD: ');
	arord=0;

	
	if ~exist('D'); % A,B need to be scaled, D,E created
		A(:,2:13)=A(:,2:13) * 0.01;  % scale ppt
		B(:,2:13)=B(:,2:13) * 0.1;   % scale temp

		% build D, E

		D= [A(3:m1,1)  A(1:m1-2,2:13)  A(2:m1-1,2:13)  A(3:m1,2:13)];
		E= [B(3:m1,1)  B(1:m1-2,2:13)  B(2:m1-1,2:13)  B(3:m1,2:13)];

	end;  % if for scaling ppt, temp and creating lagged arrays

	% Calc desired rows of D, E, C for current tree-ring series

	L1 = (D(:,1) >= a(1) & D(:,1) <= a(2)); % approp years 
	L2 = (C(:,1) >= a(1) & C(:,1) <= a(2)); % approp years 



%****************   CHECK THAT YEAR RANGE IS CONSISTENT   *****************

if (C(L2,1) ~= D(L1,1));
	error ('L2 AND L1 SPECIFY DIFFERENT YEARS IN C THAN IN D')
end

check=C(L2,1);
d1=(check(1):check(length(check)))';
if d1 ~= check 
	error ('L2 and L1 DO NOT GIVE UNBROKEN STRETCH OF YEARS')
end



%***********   WHITEN TREE-RING INDEX   ******************


C1 = C(L2,d+1);  % This is unwhitened version of tri.  C2 is
%		the whitened version.

% Model as ar(up-to-order-nhi) to get ar order.  This model will be
% built on years pointed to by L2.  The first arorder residuals will
% not be valid.  

[C2,arord,varrat,arcs]=whit1(C1,nhi,1);  % Find order for ar model.

if arord ~= 0;  %  Fit same order ar model, but extend time series
%	arord years on front end so that get valid resids for
%	full length of C2

	nhi=arord;
	j=find(L2==1);
	L3=L2;
	L3(j(1)-arord:j(1)-1) = ones(arord,1);
	[C2,arord,varrat,arcs]=whit1(C(L3,d+1),nhi,2);
	C2(1:arord)=[];	

else;  % ar order is zero -- null model
%	Whit1 will have returned arord=0, C2= C1, and the
%	(un-used) ar coefs and std devs  in arcs.
	varrat = 1.0;
end

	disp(arord);
	disp(arcs);
	disp(varrat);
	plot(C(L2,1),C1,C(L2,1),C2);
	title('ORIGINAL AND WHITENED TRI');
	pause



%****************   FORM SUB-ARRAYS OF CLIMATE SERIES  *********


% D and E hold the 37-column lagged precip and temperature series 
% covering a period larger than that for analysis.  We want a
% 12-col arrays covering only the period for analysis.

% Let DS and ES be the corresp sub-arrays for D and E.

lastmo = 37 - (12-endmo);  % tie ending month to column of D, E

DS= D(L1,lastmo-11:lastmo);
ES = E (L1,lastmo-11:lastmo);



%***********   WHITEN THE PPT AND TEMPERATURE ARRAYS  *****************

% DS and ES hold the ppt and temperature arrays for the analysis period.
% Let DW and EW hold the corresp arrays whose columns are whitened  by
% ar models of maximum possible order nhi.  Note that the first maxord
% values will be invalid for some series, where maxord is the highest 
% order model fit to any of the 24 climate series.  Need to subsequently
% truncate the front maxord years off DS, ES, DW, EW, and other arrays
% before doing correlation analysis.

maxord=0;   % initialize maximum ar order of climate variables.
parord=zeros(1,12); % Initz ar orders for ppt series.
pct1  =zeros(1,12); % Initz ar model pct var expld ppt.
tarord = zeros(1,12); % Likewise for temperature
pct2 = zeros(1,12);  % Likewise for temperature

DW=zeros(sum(L1),12);
EW=zeros(sum(L1),12);

for i = 1:12;
	[DW(:,i),parord(i),ddd,aaa]= whit1(DS(:,i),nhi,1);
	if parord(i)==0;
		pct1(i)=0;
	else
		pct1(i) = 100 * (1-ddd);
	end

	[EW(:,i),tarord(i),ddd,aaa]= whit1(ES(:,i),nhi,1);
	if tarord(i)==0;
		pct2(i)=0;
	else
		 pct2(i) = 100 * (1-ddd);
	end
	
end


%*************  TRUNCATE LEADING YEARS   ***************

maxord = max([parord tarord]);  % highest order ar model fit to ppt or 
	% temp

DS(1:maxord,:)=[];
ES(1:maxord,:)=[];
DW(1:maxord,:)=[];
EW(1:maxord,:)=[];
C1(1:maxord,:)=[];
C2(1:maxord,:)=[];

anew= [a(1)+maxord  a(2)];    % Begin and end years for correlation
ncorsize = anew(2)-anew(1) +1;  % length of time series for final correlations


%********   TWO-STANDARD ERROR BARS FOR CORRELATION  ***************

% Apply to whitened series.  Following Chatfield (1975, p. 173), the
% values falling outside the interval +- 2/sqrt(ncorsize) are signif
% at 95 % level.

sigcor= 2.0 *  1 / (sqrt(ncorsize));

%*********   CORRELATION COEFFS, BEFORE DECOUPLING  ****************

RD1=zeros(1,12);
RD2=zeros(1,12);
RE1=zeros(1,12);
RE2=zeros(1,12);

for i=1:12
	r1=corrcoef(C1,DS(:,i));
	r2=corrcoef(C2,DW(:,i));
	RD1(i) =r1(1,2);
	RD2(i) = r2(1,2);
	r1=corrcoef(C1,ES(:,i));
	r2=corrcoef(C2,EW(:,i));
	RE1(i)=r1(1,2);
	RE2(i)=r2(1,2);
end


%**********  PARTIAL CORRS YP.T AND YT.P  ******************
%
% All series whitened before analysis.
% r between tree and ppt, with effect of temp removed from both.
% r between tree and temp, with effect of ppt removed from both.
%
% Successively regress each of 24 ppt series against its paired
% monthly temperature variable.  
% in a stepwise procedure.  Let variables enter according to size of
% correlation coef of remaining variables with residuals from regression
% against all in equation so far.  Stop process when adjusted R-squared 
% is maximum.  Store the (1) number of predictors, (2) R-squared for 
% regression, and (3) Corresponding F-ratio for each regression.  Also 
% regress the tree -ring variable against the same sets of 23 potential
% predictors.  Store residuals from all regressions.

stats=zeros(1,3);
npreda=zeros(1,24);
RSQa=zeros(1,24);
FRATa=zeros(1,24);

npredat=zeros(1,24);
RSQat=zeros(1,24);
FRATat=zeros(1,24);

X=[DW EW];    % Whitened climate series of ppt, temp are potential 
		% predictors

for i = 1:24;  %  Loop for each of 24 monthly climate series.
	 if i<13
		I1=i+12;  % turn on pointer to other climate vbl same month.
	else
		I1=i-12;
	end
	y=X(:,i);    %  Predictand

	[I2,stats,c,eep(:,i),yhat] = stepr(X,y,I1);

%	I2 is pointer to I1 telling which potential predictors selected.
%  stats [1 2 3] is  R-squared, adjusted R-squared, and F-ratio
%  c are coefficients, constant first
%  eep are regression residuals
%  yhat are predicted values
	
	npreda(i) = length(I2);  % Number of predictors in final model.
	RSQa(i) = stats(1);   % R-squared
	FRATa(i) = stats(3);   % F-ratio for signif of overall equation

	y = C2;  % Now change the predictand to be the whitened tree-ring
		% index.

	[I2,stats,c,ecp(:,i),yhat] = stepr (X,y,I1);

	RSQat(i) = stats(1);
	npredat(i) = length(I2);
	FRATat(i)= stats(3);
%  ecp(:,i) are the residuals from regression of whitened tree-ring index
%  	against temp, if p is the predictor, against p, if temp is .

	disp(i)
end


%***********   CORRELATION BETWEEN REGRESSION RESIDUALS   ************
%
% Want r between pairs of regression residual series.  First member of a
% pair is residual of tree-ring index against 23 climate variables.   
% Second member of pair is residual of the left-out climate series against
% the same 23 other climate series.


RR=zeros(1,24);
for i = 1:24
	r3= corrcoef(ecp(:,i),eep(:,i));
	RR(i) = r3(1,2);
end
RD3=RR(1:12);  
RE3=RR(13:24);





%*************   DECOUPLE   BY LINEAR REGRESSION  *******************
%
% Successively regress each of 24 climate variables against the other 23
% in a stepwise procedure.  Let variables enter according to size of
% correlation coef of remaining variables with residuals from regression
% against all in equation so far.  Stop process when adjusted R-squared 
% is maximum.  Store the (1) number of predictors, (2) R-squared for 
% regression, and (3) Corresponding F-ratio for each regression.  Also 
% regress the tree-ring variable against the same sets of 23 potential
% predictors.  Store residuals from all regressions.


X=[DW EW];    % Whitened climate series of ppt, temp are potential 
		% predictors

for i = 1:24;  %  Loop for each of 24 monthly climate series.
	I1 = ones(1,24);  % pointer to climate series in X.
	I1(i) = 0;   %  Turn off pointer to the current predictand.
	y=X(:,i);    %  Predictand

	[I2,stats,c,ee(:,i),yhat] = stepr(X,y,I1);

%	I2 is pointer to I1 telling which potential predictors selected.
%  stats [1 2 3] is  R-squared, adjusted R-squared, and F-ratio
%  c are coefficients, constant first
%  ee are regression residuals
%  yhat are predicted values
	
	npred(i) = length(I2);  % Number of predictors in final model.
	RSQ(i) = stats(1);   % R-squared
	FRAT(i) = stats(3);   % F-ratio for signif of overall equation

	y = C2;  % Now change the predictand to be the whitened tree-ring
		% index.

	[I2,stats,c,ec(:,i),yhat] = stepr (X,y,I1);

%  ec(:,i) are the residuals from regression of whitened tree-ring index
%  	against all climate variables except variable i.

	disp(i)
end


%***********   CORRELATION BETWEEN REGRESSION RESIDUALS   ************
%
% Want r between pairs of regression residual series.  First member of a
% pair is residual of tree-ring index against 23 climate variables.   
% Second member of pair is residual of the left-out climate series against
% the same 23 other climate series.


RR=zeros(1,24);
for i = 1:24
	r3= corrcoef(ec(:,i),ee(:,i));
	RR(i) = r3(1,2);
end
RD4=RR(1:12);  
RE4=RR(13:24);


end;  % k2 for loop for each series 

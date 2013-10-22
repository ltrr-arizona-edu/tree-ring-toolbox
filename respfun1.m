function respfun1(P,T,X,yrs,kopt,nmos,endmo,prewh,qmax,vlab,I1,slab,Clr)
% respfun1:  response function summary of linear relationships between tree rings and climate
% CALL: respfun1(P,T,X,yrs,kopt,nmos,endmo,prewh,qmax,vlab,I1,slab,Clr);
% Last revised: 7-11-99
%
% A suite of linear analyses to identify the monthly or seasonal climate variables
% most strongly related to a tree-ring series.  Graphical summaries generated include: <BR>
% 1) correlation with monthly climate<BR>
% 2) correlation with seasonal climate<BR>
% 3) crossvalidated multiple linear regression of tree rings against seasonal climate<BR>
% 4) regression weights of tree rings on PCs of seasonal climate <BR>
% 5) regression weights from (4) in terms of original seasonal climate variables (response function)<BR>
% 6) weights of climate PCs on the seasonal climate variables<P>
%
% The degrees of freedom for CI on monthly and seasonal correlations is adjusted for autocorrelation.  
% The MLR with crossvalidation can be used to decide which climatic seasons are important to 
% the tree-ring index.  The response function is based on a principal components regression of 
% tree-ring index against PCs of seasonal climate.  Only those climate PCs whose regression coefficients 
% pass a t-test are used in computing the response function. <P>
%
% The input specifications are complicated, but there is an easy way out.  
% To avoid having to manually prepare the many input variables, you can run respfun0.m
% to automatically build a .mat file of the required input for respfun1.m.  Then simply 
% load the .mat file and type "eval(call)" at the matlab prompt to rung respfun1.m.
%
%*** IN 
%
% P(? x 13)r  monthly precipitation matrix, year in col 1
% T(? x 13)r  monthly temperature matrix, year in col 1
% X(? x 2)r   year and tree-ring index value
% kopt(1 x ?)i  Options
%   kopt(1) Accum or Average over months in forming seasons
%      ==1 accum over variable P, average over variable T
%      ==2 accum over both (Precip and Heating deg days)
%      ==3 ave over both (e.g., soil moisture and T)
%   kopt(2) Method for confidence limits on response function weights
%      ==1 using statistical distributions (t-dist)
%      ==2 bootstrap -- not yet implemented
%      ==3 index sequential, based on shifting tree-ring data relative to  climate
%          data -- not yet implemented
%      ==4 Monte Carlo: persistence AR modeled; normal assumption for residuals -- 
%          not yet implemented
%   kopt(3)  Color or B/W figures
%      ==1 color
%      ==2 B/W
%   kopt(4) Optional removal of cross-correlation between the 2 climate variables
%      ==1 Climate variables P, T  (no uncoupling)-- only available option so far
%      ==2 Climate variables P, T.P (T variable adjusted for correlation with P)
%      ==3 Climate variables P.T, T (P variable adjusted for correlation with T)
% yrs (? x 2)i  start, end year of desired analysis period for response function
%    These years are for the "tree-growth" year
% nmos (1 x 1)i number of months in the monthly window to consider for plausible
%   response (e.g., 14 if previous July through current August)
% endmo (1 x 1)i ending month (1==Jan, 12==Dec) of nmos period
% prewh (1 x ?)s  Prewhiten ('Yes') or not ('No')the tree-ring series before
%   relating to climate.  This prewhitening is done by fitting an AR(up-to-qmax)
%   model to the tree-ring data for the period ovelapping climate data
% qmax(1 x 1)i  maximum AR to consider in prewhitening tree-ring series
% vlab(2 x1)s  One-letter variable labels (e.g., ['P'; 'T'])
% I1 (? x 2)i  start and end index (relative to 1:nmos) of "seasonal"
%   groupings for the seasonal response function analysis. 
% slab (? x ?)s  ?-letter season labels.  Row size of slab will equal row
%   size of I1, or the number of seasons.  Col size of slab will be number
%   of letters in the seasonal names.  Example
%   ['Fallprev'; 'Wintcurr'; 'Summcurr'], or ['Fall'; 'Wint'; 'Summ']
% Clr -- color "triplet" matrix controlling color scheme for correl plots
%   Rows 1-3 for .05, .10 and other for color verion of P
%   Rows 4-6 for .05, ..... of T
%   Rows 7-9 for ....Black and White 
%
%
%*** OUTPUT
%
% No output arguments, just graphics
%
%
%*** REFERENCES 
%
% Tree-ring "response function" and graphical summary of monthly and seasonal
% climate signal in tree rings are discussed by Fritts, H. C. 1976. Tree Rings and 
% Climate. Academic Press. 
%
% Adjustment of degrees of freedom of correlation coefficient for autocorrelaton
% after World Meterorological Organization, 1966, Technical Note No. 79: Climatic 
% Change, WMO-No, 195.TP.100, Geneva, 80 pp.
%
% Critical 0.05 and 0.01 alpha levels for correlation coefficient after Panofsky, H.A., 
% and Brier, G. W., 1968; Some applications of statistics to meteorology, The
% Pennsylvania State University Press, p. 92.
%
% Principal components regression, including t-test for significance of regression 
% coefficients and method for back-transformation of regression weights in terms as 
% coefficients on original variables, after Mardia, K., Kent, J., and Bibby, J., 1979, 
% Multivariate Analysis: Academic Press, p. 518.
%
%*** UW FUNCTIONS CALLED (/ = subfunction)
%
% acf$$
% crospull$$
% eigen1$$
% prepmos /sideby$$
% rederr
% stepvbl1$$
% whit1$$
% /respfn1b,  /respfuna ,  /respfn1a, /subfcn01, /subfcn03 , /subfcn04
%
%*** TOOLBOXES NEEDED
%
% statistics
% system identification
%
%
%*** NOTES
%
% P and T need not be precipitation and temperature.  More generally, they can
% be any two types of monthly climatic/hydrologic variables you want to 
% relate to the tree-ring seres
%
% X:  X can be either the standard, residual or ARSTAN version of the tree-ring 
% index.  Just keep in mind that if the residual version, the AR modeling to 
% remove dependence of tree-ring on previous years' values will have little
% effect -- since this in a way has already been done
%
% yrs:  By tree-growth year, mean that the year of the tree-ring variable X.
% This year will normally be the year of the ending month of the climate data.
% For example, if period is 13-month July to July, yrs refers to the year of 
% the ending July
%
% Make sure that specified yrs is consistent with time coverage of P, T, 
% and the tree-ring data, X.  And in doing so, take into account the need for leading
% months of climate for the "first" year, and the loss of years of tree-rings
% in getting the residual series by AR modeling.  For example, say you have the 
% following data:
%
% P covers 1900-1996
% T covers 1905-1906
% X covers 1905-1985
% qmax == 2
% nmos = 18 ... 18 month period considered
% endmo = 9 .... ending month of September of growth year
%
% For this setup, yrs == [1905 1985] will not work because first avail
% residual tree-ring value for AR(2) model is in 1907.  Even if tree-ring
% data started earlier, still cannot use 1905 as start because first 
% available full climate period is 18-month period ending in 1906.  
%-------------------
%
% qmax: AR models up to order qmax will be fit for the tree-ring data for the
% full period covered by X.  The AIC is used to pick the best model. That
% model is used in prewhitening
%
%I1:  Say endmo==9 and nmos=14.  Say want 3 seasons:  Prev Fall (Aug-Nov),
%  Winter (Dec-Mar), and Spring/Summer (Apr-Sept).  Then would set
%  I1==[1 4; 5 8; 9 14]
%
% PCA regression equations coded from Mardia et al., 1979, p. 244-246
%
% Significance of correlation coefficients.  Source: Panofsky and Brier (1968, p. 92)
% Under assumption that populations are normal with zero correlation (the null hypoth),
% if N is large, the correl coef between two samples is approx normal with a mean of
% zero and a standard dev of 1/sqrt(N-2)
%
% User prompted for whether to crossvalidate MLR or skip crossvalidation
%
%------ SUBFUNCTIONS
%
% respfn1a MLR graphics
% respfn1b PCR analysis and graphics
%
% **************** OUTPUT SUMMARY ***************8
%
% FIGURE WINDOWS
%  1 monthly correlation analysis
%  2 seasonal correlation analysis 
%  3 Seasonal MLR: R-sq and crossvalidation RE vs step 
%  4 PCR -- bar chart of stdzd weights as function of step
%  5 PCR -- response function; same design as window 3
%  6 reserved for plot of selected climate PC from interactive menu
%  7 OE color map--P reserved
%  8 OE color map--T reserved
%  9...optional graphical OE summaries: e.g., B order, F order
%
% TEXT FILES
%  out?.txt (user prompted for ?)
%     Analysis period 
%     Months
%     
%     The two climate variables
%     Whether accum or average used for climate variables
%     Method(s) for condidence limits
%     Color vs B/W desired
%     Color triplet matrix
%     Whether analysis on standard or residual chron
%


kcross=questdlg('Suffer through crossvalidation of the seasonal MLR?');



%**********  BUILD MONTH LABEL COLUMN VECTOR
montxt = ['JFMAMJJASOND'];
montxt = [montxt montxt montxt];
txtsp = 24+endmo; % ending month letter
txtgo = (24+endmo) - (nmos-1);
mlab = (montxt(txtgo:txtsp))';

%****************** MAKE YEAR POINTERS TO P, T, X

% Find out if month window crosses year boundary.  Need to know because
% if does, will need monthly climate data for calendar year preceding
% yrs(1).
if nmos>endmo;
   yr1clim = yrs(1)-1;
else
   yr1clim = yrs(1);
end

% compute first year of tree-ring data needed to supply yrs(1) tree ring
% when considering need for leading values for prewhitening
yr1tree = yrs(1)-qmax;


%************ YEAR VECTORS FOR INPUT MATRICES
yrx = X(:,1);
yrP = P(:,1);
yrT = T(:,1);


%************* CHECK THAT P,T, X HAVE REQUIRED YEAR COVERAGE

if yrx(1)>yr1tree | max(yrx)<yrs(2);
   strtmp1 = sprintf('%4.0f-%4.0f',yrs);
   strtmp2 = sprintf('%4.0f-%4.0f',min(yrx),max(yrx));
   strtmp3 = sprintf('%3.0f',qmax);
   disp(['Analysis period  = ' strtmp1]);
   disp(['X coverage = ' strtmp2]);
   disp(['Max AR order for prewhitening = ' strtmp3]);
   error('Insufficient time coverage by Tree Ring matrix X (see above)');
end

if yrP(1)>yr1clim | max(yrP)<yrs(2);
   strtmp1 = sprintf('%4.0f-%4.0f',yrs);
   strtmp2 = sprintf('%4.0f-%4.0f',min(yrP),max(yrP));
   strtmp3 = sprintf('%3.0f, %3.0f',nmos, endmo);
   disp(['Analysis period  = ' strtmp1]);
   disp(['P coverage = ' strtmp2]);
   disp(['nmos, endmo = ' strtmp3]);
   error('Insufficient time coverage by climate variable P (see above)');
end

if yrT(1)>yr1clim | max(yrT)<yrs(2);
   strtmp1 = sprintf('%4.0f-%4.0f',yrs);
   strtmp2 = sprintf('%4.0f-%4.0f',min(yrT),max(yrT));
   strtmp3 = sprintf('%3.0f, %3.0f',nmos, endmo);
   disp(['Analysis period  = ' strtmp1]);
   disp(['T coverage = ' strtmp2]);
   disp(['nmos, endmo = ' strtmp3]);
   error('Insufficient time coverage by climate variable T (see above)');
end

strpd = sprintf('%4.0f-%4.0f',yrs); % string years for analysis. e.g., 1900-1980 

%****************** ORGANIZE LAGGED MONTHLY CLIMATE MATRIX
[A,mlabs]=prepmos(P,T,yrs,endmo,nmos); % Organize the monthly data,
[mA,nA]=size(A);
   
   
%************** PULL NON-WHITENED TREE-RING INDEX AND THE
%     EXTENDED INDEX NEEDED TO GET PREWHITENED INDEX
L1tree = yrx>=yrs(1) & yrx<=yrs(2);
L2tree = yrx>=yrs(1)-qmax & yrx<=yrs(2);
x1 = X(L1tree,2);   % segment of tree-ring index for climate-calibration period
xar = X(L2tree,2);  % segment of tree-ring index plus leading values needed to 
%    get AR residuals covering climate-calibration period
xarlong = X(:,2);  % index, full length


%*********  PREWHITEN TREE RINGS (OPTION) USING MODEL FIT
%      TO THE DATA ANALYSIS PERIOD COVERED BY yrs4

% Residual tree-ring series, from AR model fit to (approximately) climate calib period
[x2,k1,vrat,arcs] = whit1(xar,qmax,1); % fit ar models up to order maxar
x2(1:qmax)=[]; % AR residuals, from fit of order AR(k1) model
%      vrat is ratio of residual to original variance

% Full-length residual tree-ring series, from AR model fit to  full-length series
[x3,k13,vrat3,arcs3]=whit1(xarlong,qmax,1);
x3=x3(L1tree);

% Compute ACFs of AR-residual tree-ring indices and store results
AA = repmat(NaN,14,2);
[r2,SE22,r295]=acf(xar,5); % SE22 is two standard error bands; r292 if 95% sig level of r(1)
[r3,SE23,r395]=acf(xarlong,5);
AA(1,:)=[min(yrx)  yrs(1)];
AA(2,:)=[max(yrx)  yrs(2)];
AA(3,:)=[k13 k1];
AA(4,:)=[(1-vrat3) 1-vrat];
AA(5:9,:)=[r3' r2'];
AA(10:14,:)=[SE23' SE22'];


% Analysis on AR-filtered tree rings or not
if strcmp(prewh,'Yes');
   x1=x2; % set x1 to the prewhitened version of the calibration-period index
else
   x1=x1; % set x1 to the original (not prewhitened) index
end
  
  
%*****************  SIMPLE CORRELATION ANALYSIS, TREE VS EACH
%                   MONTHS CLIMATE

C=[x1  A(:,2:nA)]; % x1 is tree-ring index; first cols of A are monthly data for
%  the "P" variable; next cols are for the "T" variable

R1=corrcoef(C); % compute correlation mtx, tree-ring index and monthly climate variables
[mR,nR]=size(R1);
R=R1(1,2:nR);   % pull simple correl coefs, tree ring with various mos climate

% Compute 2-SE confid limits for correl coefs, adjusting for
% persistence.  
ry=zeros(1,nA-1);  % Preallocate
NN=zeros(1,nA-1);
cl=zeros(1,nA-1);
N=mA;  % Unadjusted sample size (number of years in A)
[r,SE2,r95]=acf(x1,2);
rx=r(1);   % r1 for tree-ring index
rxx=max([0 rx]); % a negative rx should not lead to any adjustment
for i=2:nA;   % Loop for each months climate
	[r,SE2,r95]=acf(A(:,i),2);
	ry(i-1)=max(abs([r(1) rx]));
	ryy=ry(i-1);
	NN(i-1)=((1-ryy*rxx) / (1+ryy*rxx)) * N;
	if NN>N, NN=N; end;  % Dont increase degrees of freedom for neg ar coef
	cl(i-1)=2.0 / sqrt(NN(i-1));
end
% Call matlab stat function norminv to get level of correl different from
% zero at the 99% and 95% levels (Panofsky and Brier , p. 92).
crit99=norminv([.995 ],0,1./sqrt(NN-2));
crit95=norminv(.975,0,1./sqrt(NN-2));
% Logical pointer to correl coefs signif at 95%, 90%, and not at 90%
Lcrit = logical(zeros(3,nmos*2));
Lcrit(3,:) = abs(R) < abs(crit95);
Lcrit(2,:) = abs(R) >= abs(crit95) & abs(R)<crit99;
Lcrit(1,:) = abs(R)>=crit99;

%--- CORRELATION PLOT
fwind=1;

% Call respfuna.m to generate plot. Note that R transposed because need col vector,
% C specifies color map, and kopt(3) controls color vs B/W.
%
respfuna(fwind,vlab,mlab,R',Lcrit,Clr,kopt(3),yrs)



%*********   SEASONAL CORRELATION ANALYSIS

nseas = size(I1,1); % number of seasons
B=repmat(NaN,mA,nseas*2); % to store seasonal data
yrA = A(:,1);

% First climate varible (P)
Atemp=A(:,2:(nmos+1)); % temporary store monthly data
for n = 1:nseas;
   ithis = I1(n,:);
   imonths=(ithis(1):ithis(2));
   num1 = length(imonths); % number of months in season
   Htemp=Atemp(:,imonths);
   % Compute col vector of sum over months
   if num1>1;
      sum1 = (sum(Htemp'))';  % cv of sum over months in season
   else; % only one month in season
      sum1=Htemp;
   end
   % Store seasonalized climate variable in B
   if kopt(1)==1 | kopt(1)==2;  % if P variable to be accumulated
      B(:,n)=sum1;
   else; % if P to be averaged over months
      B(:,n)=sum1/num1;
   end
end

% Second climate varible (T)
Atemp=A(:,(nmos+2):(2*nmos+1));
for n = 1:nseas;
   ntemp = n+nseas;
   ithis = I1(n,:);
   imonths=ithis(1):ithis(2);
   num1 = length(imonths); % number of months in season
   Htemp=Atemp(:,imonths);
   % Sum the climate variable over months
   if num1>1;
      sum1 = (sum(Htemp'))';  % cv of sum over months in season
   else;
      sum1=Htemp;
   end;
   
   % Take average rather than sum, if desired
   if kopt(1)==2;  % if T variable to be accumulated
      B(:,ntemp)=sum1;
   else; % if T to be averaged over months
      B(:,ntemp)=sum1/num1;
   end
end
% slap on year column
B=[yrA B];
[mB,nB]=size(B);

C=[x1  B(:,2:nB)]; % correlation matrix, tree ring and seasonal clim

R1=corrcoef(C);
[mR,nR]=size(R1);
R=R1(1,2:nR);   % simple correl coefs, tree ring with various seas climate

% Compute 2-SE confid limits for correl coefs, adjusting for
% persistence.  Adjust only if first order autocorr coeff positive. Base adjustment
% on the maximum r1 of the climate and tree-ring series for the pair
ry=zeros(1,nB-1);  % Preallocate
r1seas = repmat(NaN,1,nB-1);  % to hold first order autocorr of seasonal clim vbls
NN=zeros(1,nB-1); % ... for adjusted sample size
N=mB;  % Unadjusted sample size (number of years in B)
[r,SE2,r95]=acf(x1,2);
rx=r(1);   % r1 for tree-ring index
for i=2:nB;   % Loop for each season climate
   [r,SE2,r95]=acf(B(:,i),5);
   r1seas(i-1) = r(1);  
   ry(i-1)=max([ 0 r(1)]); % maximum of 0 and the first-order r of seas climate
   ryy=ry(i-1); % convenient renaming
   rxx=max([0 rx]); % effective first order autocorr of tree rings for adjustment
   numer1=(1-ryy*rxx);
   denom1=(1+ryy*rxx);
	NN(i-1)= N * numer1/denom1;
	if NN>N, NN=N; end;  % Dont increase degrees of freedom for neg ar coef
end
% Call matlab stat function norminv to get level of correl different from
% zero at the 99% and 95% levels (Panofsky and Brier , p. 92).
crit99=norminv([.995 ],0,1./sqrt(NN-2));
crit95=norminv(.975,0,1./sqrt(NN-2));
% Logical pointer to correl coefs signif at 99%, 95%, and not at 95%
Lcrit = logical(zeros(3,nseas*2));
Lcrit(3,:) = abs(R) < abs(crit95);
Lcrit(2,:) = abs(R) >= abs(crit95) & abs(R)<crit99;
Lcrit(1,:) = abs(R)>=crit99;

%--- CORRELATION PLOT
fwind=fwind+1;
% Call respfuna.m to generate plot. Note that R transposed because need col vector,
% C specifies color map, and kopt(3) controls color vs B/W.
%
respfuna(fwind,vlab,slab,R',Lcrit,Clr,kopt(3),yrs)


%******************  SEASONAL MLR
%
% 

clc;
if strcmp(kcross,'Yes'); % 
   disp(['Starting Cross-validation seasonal MLR -- and it takes time']);
end;


%-----  Build seasonal variable names (e.g., P-Fall)
nseasv = nseas*2; % number of seasonal variables
vtemp = repmat(vlab(1),nseas,1);
stemp = char(slab);
svlab = [vtemp repmat('-',nseas,1) stemp];
vtemp = repmat(vlab(2),nseas,1);
svlab=[svlab; [vtemp repmat('-',nseas,1) stemp]];
clear vtemp stemp
%  svlab is a char array of seasonal climate variable names, 1 name per row

%------ This analysis is done on standardized (zero mean, unit std dev) series
% Get the z-scorres
Z=zscore(B(:,2:nB));  % standardize lagged clim array 
[mz,nz]=size(Z);
y=zscore(x1);   % Get desired tree-ring series, and standardize

%------ Cut-off of model will be decided from cross-validation
Luse = crospull(mz,0,0); % pointer to non-left out years for crossv
% Each col of Luse defines the years to be omitted and kept for calibration
% and the year to be left out for cross-validation
% A leave-one-out strategy is used

npotent = nseasv; % number of potential predictors

%------------ Initialize storage for crossvalidation variables
RMSEv = repmat(NaN,npotent,1); % validation root-mean-square error
RMSEc = repmat(NaN,npotent,1); %  calibration root-mean-square error
RE = repmat(NaN,npotent,1); % reduction of error, from cross validation
RSQ  = repmat(NaN,npotent,1); % R-squared (calibration)
RSQp = repmat(NaN,npotent,1); % p-value of over F for equation

alpha = 0.05; % for confidence interval around the coefficients
Lin=logical(zeros(1,npotent)); % initialize pointer to potental variables in equation
Lmask = logical(zeros(1,npotent)); % do not mask any variables from consideration
iord =repmat(NaN,1,npotent); % keeps track of order of entry
nin=0; % number of predictors now in model

for n = 1:npotent; % begin cross-validation, entering 1, 2, ... npotent variables
   % Get pointer to col of Z (predictor variables) indicating the variable that
   % along with the predictors already in the equation, gives highest R-squared
   
   disp(['  ' int2str(n) ' variables in model']);
   %pause
   ipick = stepvbl1(y,Z,Lin,Lmask);
   Lin(ipick)=1;
   
   Zds =Z(:,Lin);
   [b,bint,e,eint,s] = regress(y,[ones(mz,1) Zds],alpha);
   nin=nin+1;
   iord(ipick)=nin;
   yhshort = [ones(mz,1) Zds] * b; % predicted tree-rings, calib period
   eshort = y -yhshort; % residuals
         
   % Compute standard error of estimate
   %  (See Weisberg, p. 44, 1985). Sqrt of the residual mean square)
   seest=  sqrt(sum(eshort .* eshort)/(mz-1-nin));
   
   RSQ(nin)=s(1);
   RSQp(nin)=s(3);
       
   % square root of mean square error or regression
   RMSEc(nin) =   sqrt( sum(eshort .* eshort)/(mz-1-nin));
            
   % Cross validation
   IX = crospull(mz,0,0); % pointer to crossvalidation rows of yc
   yhcv = repmat(NaN,mz,1); % to store crossvalidation predictions
   
   if strcmp(kcross,'Yes'); % if wanted crossvalidation
      for ncv = 1:mz; % loop over crossvalidation models
         ix = IX(:,ncv);
         nsum = sum(ix); % number of years of calibration data for this model
         ytemp = y(ix); % predictand
         Zdcv = [ones(nsum,1) Zds(ix,:)]; % predictor matrix
         bcv = regress(ytemp,Zdcv);
         yhcv(ncv) = [1 Zds(ncv,:)] * bcv; % predicted value for this year
      end
      ecv = y - yhcv; % crossvalidation error time series
      RMSEv(nin) =     sqrt(sum(ecv .* ecv)/mz);
      [mae,rmse,re]=rederr(mean(y),mean(y),yhcv,y);
      RE(nin)=re(1);
      rmsev=RMSEv(nin);
   end; % of if strcmp(kcross...)
   
end; % of for n=1:npotent

pmark = 0.01;  % will want to mark where MLR equation signif at better than this level

% Make the plot of change in RSQ and RE with step
fwind=3;
strout=respfn1a(svlab,iord,RSQ,RSQp,pmark,RE,fwind);



%*********  LAG-0 CROSS-CORRELATIONS BETWEEN SEASONAL P AND T


rPT =repmat(nseas,1);  % allocate to hold correlation coefs
r1P = repmat(nseas,1); % first order autocorr of P variables
r1T = repmat(nseas,1); % first order autocorr of T variables
N1 = repmat(nseas,1); % to hold adjusted sample size for each correl coef
B1 = B(:,2:(2*nseas+1));
for n = 1:nseas; % Loop over seasons
   b1 = B1(:,n); % P variable
   b2 = B1(:,n+nseas); % T variable
   rtemp = corrcoef(b1,b2);
   rPT(n) = rtemp(2);
   
   % Compute adjusted sample size (adjusted for first-order autocorr)
   [r,SE2,r95]=acf(b1,1);
   r1P(n) = r(1); 
   [r,SE2,r95]=acf(b2,1);
   r1T(n) = r(1);    
   r1big = max([r1P(n) r1T(n)]); % largest autocorrel of the 2 seasonal clim series
   if r1big<=0;  % if both series have zero or negative 1st order autocor
      N1(n) = N; % no adjustment
   else
      N1(n) = ((1-r1big) / (1+r1big)) * N;
   end
      
   
end
% Call matlab stat function norminv to get level of correl different from
% zero at the 95% and 90% levels (Panofsky and Brier , p. 92).
crit99=norminv([.995 ],0,1./sqrt(N1-2));
crit95=norminv(.975,0,1./sqrt(N1-2));
Lcrit = logical(zeros(3,nseas));
Lcrit(3,:) = abs(rPT) < abs(crit95);
Lcrit(2,:) = abs(rPT) >= abs(crit95) & abs(rPT)<crit99;
Lcrit(1,:) = abs(rPT)>=crit99;




%**********  PCA OF MONTHLY CLIMATE;  PCA REGRESSION TO PREDICT TREE

% Pack baggage
datinb{1}=y; % standardized tree-ring series
datinb{2}=Z; % standardized climate data(nseasv columns)
datinb{3}=nseasv; % number of seasonsal variables = num of seasons times 2
datinb{4}=svlab; % names of seasonalized variables
datinb{5}=strpd; % string for period of analysis
treenm='DCY Douglas-fir';
datinb{6}=treenm; % char name of tree-ring chron
kconf = kopt(2); % option for method for confidence bands 
kcolor = kopt(3); % option for color vs b/w
datinb{7} = [kconf kcolor];
fwind=fwind+1;
datinb{8}=fwind; % first figure window function respfn1b.m should use
datinb{9}=Clr;
datinb{10}=vlab;
datinb{11}=slab;


%--- call function to do work
datout=respfn1b(datinb);

% Store climate PC info needed for later plot of climate PCs
V=datout{1}; % climate PC's (seasonal), each column a PC
FV1=datout{2}; % pct variance of climate in each climate PC
FV2=datout{3}; % decimal fraction of tree-ring variance in each climate PC

% Prepare a matrix that compares tree-ring variance accounted for from stepwise
% regression on original seasonal climate variables and PCs of those
cdummy=FV2'; % decimal fraction of tree variance in each climate PC
cdummy=sort(cdummy);
cdummy=flipud(cdummy);
cdummy=cumsum(cdummy);
idummy=(1:length(cdummy))';
D = [idummy RSQ cdummy]; % Decimal fraction of variance accounted for
clear cdummy;
figure(6);
clf;
strD = sprintf('%3d  %6.2f  %6.2f\n',D')
title('Text Summary');
strmore={'Decimal fraction of Tree-ring variance accounted for in  ',...
      'stepwise regression on orignal seasonal-climate variables (A)',...
      'and on PCs of those climate variables (B)',...
      '   ',...
      '        STEP    A     B'};
text(.1,.75,strmore,'VerticalAlignment','bottom','HorizontalAlignment','left');
htextD=text(.2,.7,strD);
set(htextD,'VerticalAlignment','top');

%**************   BEGIN MENU SECTION OF PLOTS

k2wh=1; % while loop control
kboard='No';

while k2wh==1;  % Stay in plot loop
   k2men=menu('Select which to view:',...
      'Bar chart of correlations, tree rings vs monthly climate',...
      'Bar chart of correlations, tree rings vs seasonal climate',...
      'Summary of MLR of tree rings on seasonal climate variables',...
      'Summary of MLR of tree rings on PCs of seasonal climate variables',...
      'Traditional response function, on PCs of seasonal climate variables',...
      'Text summary of variance explained by MLR and PCR',....
      'Climate PCs',...
      'Toggle keyboard mode',...
      'Quit');
   if k2men==1;  % Monthly correlation analysis
      figure(1);
      if strcmp(kboard,'Yes');
         eval(['keyboard;']);
      end;
   elseif k2men==2; % seasonal correlation analysis
      figure(2);
      if strcmp(kboard,'Yes');
         eval(['keyboard;']);
      end;
   elseif k2men==3; % MLR tree rings vs seasonal climate
      figure(3);
      if strcmp(kboard,'Yes');
         eval(['keyboard;']);
      end;
   elseif k2men==4; % PCR in terms of weights on climate PCs
      figure(4);
      if strcmp(kboard,'Yes');
         eval(['keyboard;']);
      end;
   elseif k2men==5; % PCR in terms of traditional response function
      figure(5);
      if strcmp(kboard,'Yes');
         eval(['keyboard;']);
      end;
   elseif k2men==6; % text summary of variance accounted for 
      figure(6);
      if strcmp(kboard,'Yes');
         eval(['keyboard;']);
      end;
      
   elseif k2men==7; % Point & click to view climate PCs
      vlist  = cellstr(int2str((1:size(V,1))'));
      k3men=menu('Choose Climate PC to View',vlist);
      vnum=k3men;
      v=V(:,vnum);
      Fvarnc=[FV1(vnum) FV2(vnum)];
      fwind=7;
      treenm=[];
      subfcn04(fwind,vlab,slab,v,vnum,Fvarnc,kopt(3),strpd,treenm);
      if strcmp(kboard,'Yes');
         eval(['keyboard;']);
      end;
   elseif k2men==8; % toggle keyboard mode
      if strcmp(kboard,'Yes');
         kboard='No';
      else;
         kboard='Yes';
      end;
   elseif k2men==9; % quit
      k2wh=0;
   end;
end; % of k2wh


%************ SUBFUNCTIONS

function strout=respfn1a(svlab,iord,RSQ,RSQp,pmark,RE,fwind)
% respfn1a:  subfunction of respfun1.m.  Plot of RE and RSQ vs step in regression
% CALL: strout=respfn1a(svlab,iord,RSQ,RSQp,pmark,RE,fwind);
%
% Meko 10-4-98
%
%******** IN (see respfun1.m for details)
%
%  svlab -- string matrix of variable labels for season P, t
%  iord (1 x ?)i  order of entry of each seasonal variable into the regression
%  RSQ (? x 1)r  coef of multiple determination
%  RSQp (? x 1)r  p-value of overall F for equation (signif of R-squared)
%  pmark (1 x 1)r  threshold p-value for marking of "significant" equation
%  RE (? x 1)r  reductio of error statistic
%     If RE all NaN, the MLR was not crossvalidated
%  fwind -- desired plot window
%
%********* OUT
%
% strout --- string matrix that can be used as word table summarizing change
%  of R-sq and RE with step in regression
% 
%   Columns are: Step , variable name, R-sq and RE


%****************  BUILD FIGURE

figure (fwind);
nstep = length(iord); % number of steps in MLR

% Check whether MLR crossvalidated
if all(isnan(RE));
   kcross='No';
else;
   kcross='Yes';
end;


% Make the ordered set of variable labels
[y,i]=sort(iord); % sorted in ascending order of entry to equation
labord = svlab(i,:);

% Mark significant equations
Lmark = RSQp<pmark;

if any(Lmark);
   tmark = find(Lmark);
   RSQmark = RSQ(Lmark);
end

t = (1:nstep)';

if any(Lmark);
   hp1=plot(t,RSQ,'o',t,RE,'x',tmark,RSQmark,'o');
else
   hp1=plot(t,RSQ,'o',t,RE,'x');
end

% Add line between markers
line1=line(t,RSQ);
set(line1,'Color',get(hp1(1),'Color'));
line2=line(t,RE);
set(line2,'Color',get(hp1(2),'Color'));

if any(Lmark);
   set(hp1(3),'Color',get(hp1(1),'Color'),...
      'MarkerFaceColor',get(hp1(1),'Color'));
end


htxt1=text(t,RSQ+.01,labord);
set(htxt1,'Rotation',90);
txt2 = ['Accuracy'];
set(gca,...
   'XLim',[0 nstep+1],...
   'YLim',[0 max(RSQ)+0.2],...
   'XTick',[1:nstep]);

xlabel('Step');
ylabel(txt2);
title('MLR: Tree Ring Index vs Seasonal Climate');

str3=['p-value \geq ' sprintf('%5.2f',pmark)];
str4=['p-value < ' sprintf('%5.2f',pmark)];



if any(Lmark);
   if strcmp(kcross,'Yes');
      legend([hp1(3); hp1(1); hp1(2)],...
         ['\it{R}\rm^2, ' str4],...
         ['\it{R}\rm^2, ' str3],...
         'RE'); 
   else;
      legend([hp1(3); hp1(1); hp1(2)],...
         ['\it{R}\rm^2, ' str4],...
         ['\it{R}\rm^2, ' str3]); 
      
   end;
   
else
   if strcmp(kcross,'Yes');
      legend('\it{R}\rm^2','RE');
   else;
      legend('\it{R}\rm^2');
   end;
   

   txtwarn = ['No equations significant at ' num2str(pmark) ' level'];
  text(0.5,max(RSQ)+0.1,txtwarn,'FontSize',14);
   
end


%************ STRING OUTPUT
strout='Step Vbl     R-sq   RE';
strout=char(strout,blanks(5));
for n=1:nstep;
   str1 = sprintf('%2.0f  ',n);
   str2 = sprintf('%s ',labord(n,:));
   str3 = sprintf('%5.2f ',RSQ(n));
   str4 = sprintf('%5.2f ',RE(n));
   strall = [str1 str2 str3 str4];
   strout = char(strout,strall);
end




function datout=respfn1b(datin);
% respfn1b: subfunction of respfun1. PCR analysis and graphics for response function
% CALL: datout=respfn1b(datin);
%
% Meko 10-6-98
%
%******* IN
%
% datin{}
%  {1}=y; % standardized tree-ring series
%  {2}=Z; % standardized climate data(nseasv columns)
%  {3}=nseasv; % number of seasonsal variables = num of seasons times 2
%  {4}=svlab; % names of seasonalized variables
%  {5}=strpd; % string for period of analysis
%  {6}=treenm; % char name of tree-ring chron
%  {7} (1 x 2)i  options
%     1=kconf  method for confidence bands around PCR weights
%         1==statistical (t-dist)
%         2==bootstrap
%         3==circular shifted
%         4==monte carlo
%  {8}=fwind  first figure window available to this function
%  {9} Clr (3 x 9)r color matrix (see respfun1.m)
%  {10} dlab (2 x 1)s data-type labels (e.g., ['P';'T']
%  {11} slab (1 x ?)c  cell vector of season names
%
%
%********* OUT
%
% datout{}
%  {1} V eigenvectors of the climate variables, each col an eigenvector
%  {2} teach (? x 1)r percentage of climate variance for each PC
%  {3} V1a (? x 1)r decimal fraction of tree variance described by each
%      climate PC
%
%
%**********  NOTES
%
%
% Effective sample size is not adjusted for autocorrelation in the tree-ring
% or climate variables.  Keep that in mind when comparing the significance
% of weights on monthly climate by this method and by the correlation method
% in respfun1.  The correlation method does optionally adjust the effective
% sample size in computation of signifance.
%
%
% Three figure windows are produced.
% First window:
%  top plot: PCR weights in bar chart, with 99% and 95% signif color marked.
%   Ordering along x axis is by number of climate PC.  
% Second Window: R-sqd and RE plot vs step, with RE based on cross-validation
%     x-axis ordered as climate PC#. Bars annotated with entrystep 
% Third Window:
%  Response function weights on monthly climate variables, from PCR model selected
%     interactively by user after he looks at first window. Third window not produced
%     until user clicks on PC in plot in second window
%
% User actions.  User looks at plots in second window and decides a cutoff number
% of climate PCs to include as predictors in the model. User allowed only
% to use PCs whose regression coeffs signif at .05 or better. User can use fewer
% PC's. In limit, might just use first PC to enter.  User tells function the last
% PC he wants to include by pointing to it's bar on the top plot in window 1.



%------- UNLOAD
y = datin{1}; 
Z = datin{2};
nseasv = datin{3};
svlab=datin{4};
strpd=datin{5};
treenm=datin{6};
kconf=datin{7}(1);
kcolor=datin{7}(2);
fwind = datin{8};
Clr=datin{9};
dlab=datin{10};
slab=datin{11};

%----- Size
[mz,nz]=size(Z);

nullfg2=0;  % this will change to 1 if no significant predictors
[RR,V,L,S,F]=eigen1(Z,1);  % pca on correl mtx of Z
W=Z*V;  %  Amplitudes of PCs of Z (cols).  Recall that Z contains standardized
% variables, and so they have mean zero.  As a result, means of cols of W
% are also zero

% Compute variance explained of climate by each climate PC

ts=sum(diag(L));
teach=100 * diag(L) ./ ts;


a=W\y;  % Estim coefs of model to predict tree from clim eig amplitudes
ep = y - W*a;  % Compute regression residuals

% Compute pct variance of tree growth explained by each climate 
% eigenvector.

Rtemp=corrcoef([y W]);
V1a=(Rtemp(1,2:nz+1)) .^2;  % proportion variance expld
V1b= sum(V1a);  % Total proportion variance expld by model with all pc's

epsq=ep' * ep;  % sum of squares of residuals of regr of tree variable on the
			  % full set of  pc amps of climate array.

term1=epsq(ones(nz,1),:); % Terms needed for eq 8.85, p. 245 in 
aa=(mz-nz-1);             % Mardia et al. 1979
term2=aa(ones(nz,1),:);
term3= sqrt(mz * diag(L));

% Find .975 and .995 probability points of "t" distribution:  for 95%
% and 99% confidence limits around regression coefficients a.
% Only the 95% cl is used in the plots, but the 99% is there if you need
% it.
t1=tinv(.975,aa);
t2=tinv(.995,aa);
t95=t1 * sqrt(term1 ./ term2) ./ term3;
%((t1 * term1 ./ term2) .^2) ./  term3;

t99=t2 * sqrt(term1 ./ term2) ./ term3;
%t99=((t2 * term1 ./ term2) .^2) ./  term3;

% Compute standardized PCR weights, corresponding to terms in eqn 8.8.5, p 245
% These will be used in plotting
aplot = (a .* term3) ./ sqrt(term1 ./ term2);


% Logical pointer to PCR stdzd regression weights signif at 99%, 95%, and not at 95%
Lc = logical(zeros(nseasv,3));
Lc(:,3) = abs(a) < abs(t95);
Lc(:,2) = abs(a) >= abs(t95) & abs(a)<t99;
Lc(:,1) = abs(a)>=t99;



J = abs(a) < abs(t95);  % Ones point to elements to zero out
Jr = ~J;  % 0/1 cv, 1s point to nonzero elements of a
ev = epsq / mz;  %  mean sum of squares of residuals from PC regression using all
%		PCs as predictors.
evstar = ev * mz/aa; % estimated variance of residuals; takes into account
%   number of df in computing variance from observed residuals
vexp=sum(V1a(Jr));  % Total proport. variance tree growth explained by the 
%	restricted set of climate eigenvectors
Jrf = find(Jr);  % Subscripts of cols of A corresp to restricted set
%      of PC predictors

%******  F-level and significance for PC regression using only
%        coefs sig at 95% level



if sum(Jr)==0, nullfg2=1; end;  % No significant PC coefficients--null model

% F-test for significance of the multiple correlation coefficient
% See Panofsky and Brier, p. 113.  In table 32 of that reference, see that
% the F is computed from R-squared, N and p.  R-squared is the proportion of
% variance accounted for by regression, which in my code is vexp.  N is the number
% of observations, which is mz; and p is the number of predictors, which is 
% sum(Jr) in my notation.


Lcrit = logical(zeros(3,nseasv));

if nullfg2==0;    %  At least one signif pc coefficient; proceed
   dfreg=sum(Jr); % degrees of freedom for regression sum of squares
   dfres=mz-sum(Jr)-1;  % deg freedom for residual sum of squares
   F2=vexp*dfres/((1-vexp)*dfreg); % Compute F ratio
   ff2 = finv(0.95,dfreg,dfres); % table value of F
   
      
   %***********  TRANSFORM THE PC-REGRESSION WEIGHTS BACK INTO 
   %***********  WEIGHTS ON THE INDIVIDUAL MONTHLY CLIMATE VARIABLES, AND
   %***********  COMPUTE ASSOCIATED ERROR BARS (MARDIA ET AL. 1979, P. 246)
   
   % Retain only those components whose  tree vs pc  regression coefs were
   % significant at 95% level.  Form a restricted least squares estimator
   % by setting all non-significant a's to zero.
   
   arest=a;  % initialize restricted estimator
   arest(J)=zeros(sum(J),1);
   b=V * arest;  %   Transformed coefs -- on the original monthly climate
   
   
   bv=zeros(length(Jr),1);
   nkeep=sum(Jr);  % number of pc s retained
   %bv=zeros(Jr);  % Preallocate
   
   Ld=diag(L);
   
   for i=1:nz;  % for each monthly climate variable
      bv(i) = (evstar / mz) * sum ((V(i,Jr) .^2) ./ (Ld(Jr))');
   end
   
   % The estimated standard error
   bvse=sqrt(bv);  % standard errors for coefs of transformed
   %		model.
   
   % Can use a t-distribution (see Draper & Smith 1981, p. 25)to test
   % significance.  A b(i) divided by its standard error is distributed
   % at t, with df the same as used to compute the error variance, which
   % is n-p-1, or aa
   
   % Call matlab stat function norminv to get level of coefs different from
   % zero at the 99% and 95% levels Mardia, p. 246
   %crit99=norminv([.995 ],0,bvse');
   crit99=tinv(.995,aa);
   %crit95=norminv(.975,0,bvse');
   crit95=tinv(.975,aa);
   % Logical pointer to correl coefs signif at 99%, 95%, and not at 95%
   Lcrit(3,:) = abs((b ./ bvse)') < abs(crit95);
   Lcrit(2,:) = abs((b ./ bvse)') >= abs(crit95) & abs((b ./ bvse)')<crit99;
   Lcrit(1,:) = abs((b ./ bvse)')>=crit99;
   
      
else;  % No significant pc-regression coefs -- null model
   clc, home;
   disp('NO SIGNIFICANT PC-REGRESSION COEFS IN MODEL');
	disp('ALL PCs TOGETHER EXPLAIN FOLLOWING PROPORTION TREE-RING VARIANCE: ');
   disp(V1b);
   disp('Press any key tO continue');
   pause;
   
end;  % of code dealing with Null model

%----- GATHER OUTPUT ARGUMENT DATA
datout{1}=V;
datout{2}=teach;
datout{3}=V1a;

%************** GRAPHICS CALL -- BAR GRAPH OF STANDARDIZED PCR WEIGHTS 
subfcn01(fwind,aplot',Lc',Clr,kcolor);

%****** GRAPHICS CALL-- BAR GRAPHS OF WEIGHTS ON ORIGINAL SEASONAL CLIMATE VARIABLES
fwind=fwind+1;

if nullfg2==0; % if at least one climate PC entered the PCR model
   subfcn03(fwind,dlab,slab,b',Lcrit,Clr,kcolor,strpd,treenm);
else; % No regression coefs significant in PCR model
   figure (fwind);
   title('NO RESPONSE FUNCTION JUSTIFIED: NO SIGNIFICANT PC PREDICTORS');
end




%********* SUBFUNCTIONS

function subfcn01(fwind,wgt,L,C,kopt)
%
% Sub-function of respfn1b.  PCR plots of stdzd weights
%
%************* IN
%
% fwind (1 x 1)i figure window
% wgt (1 x 2*nwgt)r  weights or correlations for the months or seas, P then T 
% L (3 x 2*nwgt)L  indicates whether weights signifant at .05, .10
%     row 1 applies to .01 level.  A "1" means signif, a "0" not
%     row 2 applies to .05 level. ....
%     row 3 applies to "insignificant" correls
% C (9 x 3)r  color schemes
%   Rows 1-3 for .01, .05 and other for color verion of P
%   Rows 4-6 for .01, ..... of T
%   Rows 7-9 for ....Black and White 
% kopt -----options
%   kopt(1)  color or b/w
%      ==1 color
%      ==2 b/w

figure(fwind)

statusc = close(fwind);
if statusc~=1;
   close(fwind);
end

figure(fwind);

set(gcf,'RendererMode','manual');
set(gcf,'Renderer','zbuffer');


set(gcf,'DefaultLineLineWidth',2.0);
set(gcf,'DefaultTextFontWeight','bold');

% Calculate number of weights
nwgt = length(wgt);
LP=L;

h1 = bar(wgt);
x1 = get(h1,'Xdata');
y1 = get(h1,'Ydata');


%Color patches
if kopt(1)==1; % color
   c1 =C(1,:); c2 = C(2,:); c3=C(3,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CP = zeros(1,nwgt,3);
n99 = sum(LP(1,:));
if n99>0;
   CP(1,LP(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LP(2,:));
if n95>0;
   CP(1,LP(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LP(3,:));
if nother>0;
   CP(1,LP(3,:),:) = repmat(c3,nother,1);
end
hpatch1 = patch(x1,y1,CP);
set(hpatch1,'LineWidth',2);

% Adjust parent of patch
h1 = get(hpatch1,'parent');

set(h1,'XLim',[0 nwgt+2],...
   'XTick',[1:nwgt],...
   'Xgrid','on');

% Set + y limit at twice largest absolute value of any weight -- room for legend
maxwgt=max(wgt);
yylim = get(gca,'Ylim');
set(gca,'Ylim',[yylim(1)  2*maxwgt]);
yspace = maxwgt;

xlabel('Climate PC Number');


% Horiz Line at zero weight
hold on;
line([0 nwgt+1],[0 0],'Color',[0 0 0]);
hold off;


% Reference points for positioning legend
xxlim = get(gca,'XLim');
yylim = get(gca,'Ylim');
xwide = diff(xxlim);
ywide = diff(yylim);


%**************Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
ygo1 = yylim(2)-yspace/20;
ydel = yspace/5;
ygo2 = ygo1 - 2*ydel;

% Upper part of legend
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [ygo1-ydel ygo1 ygo1 ygo1-ydel];
hlegp1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,ygo1-ydel/2,'\alpha=0.01','FontWeight','bold');

% Lower part of legend
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [ygo2-ydel ygo2 ygo2 ygo2-ydel];
hlegp2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,ygo2-ydel/2,'\alpha=0.05');

ylabel('Standardized Coefficient');
title('PCR: Weights on PC''s of Seasonal Climate Variables',...
   'FontSize',12,'FontWeight','bold');



function subfcn03(fwind,dlab,slab,wgt,L,C,kopt,strpd,treenm)
%
% Sub-function for respfun1. PCR graphs of re-transformed weights on original
% seasonal climate variables
%
%************* IN
%
% fwind (1 x 1)i figure window
% dlab (2 x 1)s data labels (e.g., ['P';'T']
% slab (1 x ?) cell vector of season names
% wgt (1 x 2*nmos)r  weights on seasonalized P then T 
% L (3 x 2*nmos)L  indicates whether weights signifant at .05, .10
%     row 1 applies to .01 level.  A "1" means signif, a "0" not
%     row 2 applies to .05 level. ....
%     row 3 applies to "insignificant" correls
% C (9 x 3)r  color schemes
%   Rows 1-3 for .01, .05 and other for color verion of P
%   Rows 4-6 for .01, ..... of T
%   Rows 7-9 for ....Black and White 
% kopt -----options
%   kopt(1)  color or b/w
%      ==1 color
%      ==2 b/w
% strpd (1 x ?)s labeling string with period of analysis
% treenm (1 x ?s  label string of tree-ring series name
%
%*********** NOTES 
%
% 


figure(fwind)

statusc = close(fwind);
if statusc~=1;
   close(fwind);
end

figure(fwind);

set(gcf,'RendererMode','manual');
set(gcf,'Renderer','zbuffer');


set(gcf,'DefaultLineLineWidth',2.0);
set(gcf,'DefaultTextFontWeight','bold');

%Compute number of seasons
nmos = size(slab,2);


%----------  P axes (at bottom)
ax1 = axes;
set(ax1,'Position',[.1 .09 .8 .4]);

% pull P weights
pwgt = wgt (1:nmos);
LP = L(:,1:nmos);


h1 = bar(pwgt);
x1 = get(h1,'Xdata');
y1 = get(h1,'Ydata');


%Color patches
if kopt(1)==1; % color
   c1 =C(1,:); c2 = C(2,:); c3=C(3,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CP = zeros(1,nmos,3);
n99 = sum(LP(1,:));
if n99>0;
   CP(1,LP(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LP(2,:));
if n95>0;
   CP(1,LP(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LP(3,:));
if nother>0;
   CP(1,LP(3,:),:) = repmat(c3,nother,1);
end
hpatch1 = patch(x1,y1,CP);
set(hpatch1,'LineWidth',2);

% Adjust parent of patch
h1 = get(hpatch1,'parent');

set(h1,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',slab,...
   'Xgrid','on');
%,...
  % 'YLim',[-1 1]);

xlabel('Season','FontWeight','bold');

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,0.70,[dlab(1,:)],'FontSize',14);
set(gca,'FontWeight','bold');

hold on;
line([0 nmos+1],[0 0],'Color',[0 0 0]);
hold off;

%Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [.75 .95 .95 .75];
hlegp1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,.85,'\alpha=0.01','FontWeight','bold');
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [.50 .70 .70 .50];
hlegp2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,.60,'\alpha=0.05');

%************** T axes (one up from bottom)
ax2 = axes;
set(ax2,'Position',[.1 .52 .8 .4]);

% pull T weights
twgt = wgt ((nmos+1):(2*nmos));
LT = L(:,(nmos+1):(2*nmos));

hb = bar(twgt);
x2 = get(hb,'Xdata');
y2 = get(hb,'Ydata');

%Color patches
if kopt(1)==1; % color
   c1 =C(4,:); c2 = C(5,:); c3=C(6,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CT = zeros(1,nmos,3);
n99 = sum(LT(1,:));
if n99>0;
   CT(1,LT(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LT(2,:));
if n95>0;
   CT(1,LT(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LT(3,:));
if nother>0;
   CT(1,LT(3,:),:) = repmat(c3,nother,1);
end

% Clear bar chart and add patch
cla;
hpatch2 = patch(x2,y2,CT);
set(hpatch2,'LineWidth',2);

% Adjust parent of patch
h2 = get(hpatch2,'parent');

set(h2,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',' ',...
   'Xgrid','on');

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,0.75,[dlab(2,:)],'FontSize',14);


hold on;
line([0 nmos+1],[0 0],'color',[0 0 0]);
hold off;

set(gca,'FontWeight','bold');

%Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [.75 .95 .95 .75];
hlegt1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,.85,'\alpha=0.01');
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [.50 .70 .70 .50];
hlegt2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,.60,'\alpha=0.05');

txt1 = [treenm ';  ' strpd ';  PCR Weights on Seasonal Climate'];
title('PCR: Weights on Original Seasonal Climate Variables',...
   'FontSize',12,'FontWeight','bold');

function respfuna(fwind,dlab,mlab,wgt,L,C,kopt,yrs)
%
% Sub-function for respfun1 that makes the page of plots
%
%************* IN
%
% fwind (1 x 1)i figure window
% dlab(2 x ?)s  label for data types (e. g., ['P';'T'])
% mlab (1 x nmos)s month labels, assumed to be in correct order; might
%   also be cell rv of season labels
% wgt (1 x 2*nmos)r  weights or correlations for the months or seas, P then T 
% L (3 x 2*nmos)L  indicates whether weights signifant at .05, .10
%     row 1 applies to .01 level.  A "1" means signif, a "0" not
%     row 2 applies to .05 level. ....
%     row 3 applies to "insignificant" correls
% C (9 x 3)r  color schemes
%   Rows 1-3 for .01, .05 and other for color verion of P
%   Rows 4-6 for .01, ..... of T
%   Rows 7-9 for ....Black and White 
% kopt -----options
%   kopt(1)  color or b/w
%      ==1 color
%      ==2 b/w
% yrs (1 x 2)i start, end year of analysis period (unadjusted N for correlations)
%
%*********** NOTES 
%
% A weights will be zero for beyond model order.  For example, if order of model
%   is 3, rows 8-9 will be zero  for the appropriate col of A
% 

% 


figure(fwind)

statusc = close(fwind);
if statusc~=1;
   close(fwind);
end

figure(fwind);

set(gcf,'RendererMode','manual');
set(gcf,'Renderer','zbuffer');


set(gcf,'DefaultLineLineWidth',2.0);
set(gcf,'DefaultTextFontWeight','bold');

if ischar(mlab);
   kmode1 = 1; % monthly model
   nmos = size(mlab,1);  % number of months in window
elseif iscell(mlab);
   kmode1 = 2; % seasonal mode
   nmos = size(mlab,2); % here mlab is a 1-row cell of season names
else;
   error('unknown mode');
end


%----------  P axes (at bottom)
ax1 = axes;
set(ax1,'Position',[.1 .09 .8 .4]);

% pull P weights
pwgt = wgt (1:nmos);
LP = L(:,1:nmos);


h1 = bar(pwgt);
x1 = get(h1,'Xdata');
y1 = get(h1,'Ydata');


%Color patches
if kopt(1)==1; % color
   c1 =C(1,:); c2 = C(2,:); c3=C(3,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CP = zeros(1,nmos,3);
n99 = sum(LP(1,:));
if n99>0;
   CP(1,LP(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LP(2,:));
if n95>0;
   CP(1,LP(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LP(3,:));
if nother>0;
   CP(1,LP(3,:),:) = repmat(c3,nother,1);
end
hpatch1 = patch(x1,y1,CP);
set(hpatch1,'LineWidth',2);

% Adjust parent of patch
h1 = get(hpatch1,'parent');

set(h1,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',mlab,...
   'Xgrid','on',...
   'YLim',[-1 1]);

if kmode1==1;
   xlabel('Month','FontWeight','bold');
else;
   xlabel('Season','FontWeight','bold');
end

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,0.70,[dlab(1,:)],'FontSize',14);
set(gca,'FontWeight','bold');

hold on;
line([0 nmos+1],[0 0],'Color',[0 0 0]);
hold off;

%Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [.75 .95 .95 .75];
hlegp1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,.85,'\alpha=0.01','FontWeight','bold');
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [.50 .70 .70 .50];
hlegp2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,.60,'\alpha=0.05');

%************** T axes (one up from bottom)
ax2 = axes;
set(ax2,'Position',[.1 .52 .8 .4]);

% pull T weights
twgt = wgt ((nmos+1):(2*nmos));
LT = L(:,(nmos+1):(2*nmos));

hb = bar(twgt);
x2 = get(hb,'Xdata');
y2 = get(hb,'Ydata');

%Color patches
if kopt(1)==1; % color
   c1 =C(4,:); c2 = C(5,:); c3=C(6,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CT = zeros(1,nmos,3);
n99 = sum(LT(1,:));
if n99>0;
   CT(1,LT(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LT(2,:));
if n95>0;
   CT(1,LT(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LT(3,:));
if nother>0;
   CT(1,LT(3,:),:) = repmat(c3,nother,1);
end

% Clear bar chart and add patch
cla;
hpatch2 = patch(x2,y2,CT);
set(hpatch2,'LineWidth',2);

% Adjust parent of patch
h2 = get(hpatch2,'parent');

set(h2,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',' ',...
   'Xgrid','on','YLim',[-1 1]);

if kmode1==1;
   %xlabel('Month');
else;
   %xlabel('Season');
end

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,0.75,[dlab(2,:)],'FontSize',14);


hold on;
line([0 nmos+1],[0 0],'color',[0 0 0]);
hold off;

set(gca,'FontWeight','bold');

%Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [.75 .95 .95 .75];
hlegt1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,.85,'\alpha=0.01');
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [.50 .70 .70 .50];
hlegt2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,.60,'\alpha=0.05');

if kmode1==1;
   txt1 = 'Monthly Correlation Analysis, ';
else;
   txt1 = 'Seasonal Correlation Analysis, ';
end
txt2 = [int2str(yrs(1)) '-' int2str(yrs(2))]; 
txt3=[txt1 txt2];
title(txt3,'FontSize',14,'FontWeight','bold');


function subfcn04(fwind,dlab,slab,v,vnum,Fvarnc,kopt,strpd,treenm)
%
% Sub-function for respfun1. Bar graph of weights of a climate PC on original
% seasonal climate variables
%
%************* IN
%
% fwind (1 x 1)i figure window
% dlab (2 x 1)s data labels (e.g., ['P';'T']
% slab (1 x ?) cell vector of season names
% v (2*nseas x 1)r seasonal-climate PC, each value a weight on a seasonal climate variable
%   Convention followed is "P" variable first, then "T" variable
% vnum (1 x 1)i  which climate pc is v?
% Fvarnc (1 x 2)r fractional percentage of climate variance accounted for
%   by this PC (Fvarnc(1)), and fractional percentage of tree-ring variable
%   accounted for in PCR by this climate PC (Fvarnc(2))
% kopt -----options
%   kopt(1)  color or b/w
%      ==1 color
%      ==2 b/w
% strpd (1 x ?)s labeling string with period of analysis
% treenm (1 x ?)s  tree-ring series name
% 
%
%*********** NOTES 

figure(fwind)

statusc = close(fwind);
if statusc~=1;
   close(fwind);
end

figure(fwind);

set(gcf,'RendererMode','manual');
set(gcf,'Renderer','zbuffer');


set(gcf,'DefaultLineLineWidth',2.0);
set(gcf,'DefaultTextFontWeight','bold');

%Compute number of seasons
nmos = size(slab,2);


% Compute y-axis limits
yeps = abs((max(v) - min(v))/20);
yaxlim=[min(v)-yeps max(v)+yeps];


%----------  P axes (at bottom)
ax1 = axes;
set(ax1,'Position',[.1 .09 .8 .35]);

% pull P weights
pwgt = v (1:nmos);

% Draw bar for "P" variable weights
h1 = bar(pwgt);
x1 = get(h1,'Xdata');
y1 = get(h1,'Ydata');

hpatch1 = patch(x1,y1,[.5 .5 .5]);
set(hpatch1,'LineWidth',2);

% Adjust parent of patch
h1 = get(hpatch1,'parent');

set(h1,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',slab,...
   'Xgrid','on');
%,...
  % 'YLim',[-1 1]);

xlabel('Season','FontWeight','bold');

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,max(v)-2*yeps,[dlab(1,:)],'FontSize',14);
set(gca,'FontWeight','bold','Ylim',yaxlim);

hold on;
line([0 nmos+1],[0 0],'Color',[0 0 0]);
hold off;

%************** T axes (one up from bottom)
ax2 = axes;
set(ax2,'Position',[.1 .47 .8 .35]);

% pull T weights
twgt = v ((nmos+1):(2*nmos));

hb = bar(twgt);
x2 = get(hb,'Xdata');
y2 = get(hb,'Ydata');


% Clear bar chart and add patch
cla;
hpatch2 = patch(x2,y2,[.5 .5 .5]);
set(hpatch2,'LineWidth',2);

% Adjust parent of patch
h2 = get(hpatch2,'parent');

set(h2,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',' ',...
   'Xgrid','on');

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,max(v)-2*yeps,[dlab(2,:)],'FontSize',14);


hold on;
line([0 nmos+1],[0 0],'color',[0 0 0]);
hold off;

set(gca,'FontWeight','bold','Ylim',yaxlim);

str1=['Weights of Climate PC #' int2str(vnum)];
str2=sprintf('  %4.1f percent of climate variance for %s',Fvarnc(1),strpd);
str3=sprintf('  %4.1f percent of Tree-Ring variance in PCR model',Fvarnc(2)*100);
str4={str1,str2,str3};

title(str4,'FontSize',10,'FontWeight','bold');





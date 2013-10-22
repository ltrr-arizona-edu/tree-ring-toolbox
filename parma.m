% arma.m

% a 386-MATLAB  program for low-order autoregressive-moving average
% modeling of time series.

% by David M. Meko, Laboratory of Tree-Ring 
% Research, University of Arizona, 85721
% August, 1990

%################  DESCRIPTION ######################################
%
% This program performs exploratory time-series modeling, considering
% the possible model as one of four candidate low-order ARMA models:
% AR(1), AR(2), MA(1), and ARMA(1,1).
%
% References are the MATLAB  manual, the book "System Identification:
% Theory for the User", by Lennart Ljung, 1987, Prentice Hall;  and
% the book "Applied Modeling of Hydrologic Time Series", by  J.D. Salas,
% J. W. Delleur, V. Yevjevich, and W.L. Lane, 1980, Water Resources
% Publications.
%
% The main steps followed by the program:
%
% 1. Read in the time series and subtract its mean.
% 2. Loop through each of the four candidate models, providing
%    graphic output to the screen to aid in model selection.
%    Compute the Akaike Information Criterion, Porte Manteau
%    statistic, and other quantities. Estimate model parameters
%    and their standard errors.
% 3. Store summary information on the selected model for possible
%    later use.  Later use can include input to a program to 
%    "whiten" series using the selected models.
% 4. Repeat  items 1-3 for each time series
%
% ######  PRELIMINARY SETUP ##########################################
%
% 1. Store in a ".mat" file three arrays:
%	1. A - the data array, with each row a year, each col a time series.
%		   Column 1 should be the year.
% 	2. colwhich - a row vector telling which time series from A
%		   are to be analyzed.  
% 	3. ends - ?x3 array, beginning and ending year on which model is 
%		   to be based for each of ? series.
%
% 2.  Start Matlab.
% 3.  >load "?.mat", whatever you have named the ".mat" file holding
%		A, colwhich and ends.
% 4.  arma.... to start the program.
%
%##################  FEATURES #######################################
%
%	*	Screen plots:
%			1. ACF of time series, with large-lag 2-se bars
%			2. AIC and Portmanteau statistic versus candidate model.
%			3. Estimated parameters with 2-standard-error bars
%			4. Time series plot of original and "whitened" data for
%			   candidate model--four frames on one screen.
%			5. Autocorrelations of model residuals
%
%	*	Optional repeated recall of various plots until decision
% 		  	is made on model selection.
%
%	*	Automatic storage of key info on selected model for each
%			series for later use.
%
%	*	Easily adapted (by inserting "keyboard" command after 
%			relevant "Plot" commands) to produce optional
%			hard-copy plots of any screen plot on LaserJet,
%			Epsom printer, or HP Plottter
%
%	*	Convenient for large-volume data analysis. Originally tested
%			on time series array 1700-1979 with 249 data series
%			or columns.  Use of "colwhich" vector makes analyzing
%			a specific subset of the variables easy.  Time period
%			of modeling for each time series can also be controlled.
%
%#################### LIST OF VARIABLES #############################
%
% acf... autocorrelation function, lags 0-24, of tme series z
% akaik ...4x1 vector of AIC values for each model tried
%
% ccf...cross correl fctn, lags 0-24, of time series z
% colwhich...1x?,  which series in A are to be analyzed.  For example
%	[5 7 9 15] designates series 5,7,9, and 15,  which would be
%	columns 6,8,10 and 16 of A.
% cume...24x1 vector used to accumulate squared acf values in computation
%	of standard error bars for acf.
% e...	residuals of an estimated AR model 
% ends... ncols X 3, tells period on which model for each series
% 	will be based.  Each row for a series.  First col is the series # 
%   or the col the series is in A, , second col 
%   is beginning valid year of that series,   and third col 
%   is  the ending valid year for modeling the series.
% evar... 4X1 vector of variance of residuals for each model.  Computed over only
%    those observations beyond the first nparam
% E...  matrix of residuals for a given time series. Each of four
%  cols corresps to one of the four candidate models.  Col 1 is the year.
% k...loop index for which series in A to analyze next.  Note that
%	k=1 corresponds to the second row of A, since row 1 is the year.
% k1..  loop index for trying different models
% k2...menu option for which plot or table to view on screen
% lag... 	[0:1:24]', vector for abscissa of plot of acf of resids 
% m...  used in various places to hold row size of an array
% m1... holds row size of ends after call to size.m
% mnp...tells frame setup for subplot.  Used in various places in pgm.
% modnum...[1 2 3 4]', id #s for each model to try
% mods... ? X 3 array of selected model order for each series
%   col 1 - series #, col2 - ar order,  col 3 - ma order
% modstat1...? X 6   summary of selected models.
%   col 1 - series #
%   col 2 - Porte Manteau Statistics based on lags 1-24
%   col 3 - Ratio of error variance to original variance
%   col 4,5,6 parm ests: ar-1 ar-2 ma-1;  zero if none
% mz... mean value of original time series z
% n...  used in various places to hold col size of an array
% n1... holds col size of ends after call to size
% ncols... number of columns in A
% ngo...subscr of yr1 telling first year to plot out in dual
% 	time series plot of original and whitened data.  
% nmods... number of different models to try
% nn... order of current model being estimated
% NN ...  orders of ARMA models to be tried .
%	Col 1 gives AR order, col 2 MA order
% nparam... totl # of arma params in a model
% nrows... number of rows in A, which equals number of years in A
%	program assumes A has been trimmed so that the rows contain the
%	valid analysis period.
% nyears ... number of valid years to analyze in series
% opt1... screen prompted control over period to be analyzed
%	Y- same period (all rows) for all  series
%	N- ends holds beginning and ending years
% parest ... 1X8, estimates of parameters of the four models tried
%	Elements 1,8 are dummy 0.   2 holds the AR(1) coef, 3 and 4 hold
%	the two ar(2) coefs, 5 holds the MA(1) coef, and 6 and 7 hold
%     the ar and ma coefs of the ARMA(1,1) model.  Taken from
%	row 3 of the th array.
% pct... 4X1 vector of ratio of residual to orig variance for each model
% q... portmanteau statistic.  See O.D. Anderson, 1976, eqn 9.10, p. 84.
%	Computed on lags 1-24.  Df for test is therefore 24 - AR order.
% Q ... vector of Porte Manteau Q for each of four models
% q23df...chi-sqd for 95% prob, 23 df.  Used in port manteau
% q22df...chi-sqd for 95% prob, 22 df.  Used in port manteau
% r... acvf of residuals of model fit for a specific series and model, as
%	returned by function resid.
% r1...acf of  resids of model fit for a specific series and model.
%	Contains acf for lags 1-24.  Computed by dividing r(2:25) by r(1).
% R...25 x nmods, autocorrelations of residuals from each model fit
%	for a given time series.
% rowgo...row subscript in A of first valid year of the "model period" of 
%  	a time series. May vary from series to series, depending on "ends" array.
%	Ends array holds actual years, while rowgo is a subscript rlative
%     to row 1 of A.
% rowstop... row subscript in A of last valid...
% seax... 1X8 , axis for plot of standard errors, ses. First 
%  	and last elements arezero, and, next elements are 1,2,3,4,5,6,7.
% seaxzero...1x8 vector of zeros for horiz line on plot of param ests
% sebar... 24x1 vector of 2-stand-errors (large-lag) for acf of time series z.
%	Corresp to error bar values for acf at lags 1-24 years  
%	Used for error bars on plot of acf.
% ses... 1X8, 2-standard errors for each parameter estimate from 
%	each model tried.  Taken from rows 4, or 4 and 5 of the theta array
%	First and last elements are dummy zeros.
% sumr2...24x1 vector of squared acf estimates for lags 1-24 years.
% TBL1...matrix model results: series #, order model, %var removed, 
%	portmanteau q
% th....packed matrix of parameter estimation results:
%  row 1 var(e), sampling interval T, nu,na,nb,nc,nd,nf,nk
%  row 2 FPE, year, mo, date, hr, min, and command model gen by
%  row 3 vector of est params: A,B,C,D, and F
%  row 4 to 3 +n:  estimated param covariance mtx, where n is sum of
%		all orders (number of est params)
% V ... loss functions (top row) for each candidate model
%       model order (row 2)
%       sample size of resids for largest order model (last col)
% varz... variance of z, computed only after deleting number of
%   initial years equal to selected AR order.  This deletion to give
%   better comparison with variance of residuals.
% VMOD... AIC for models tested
% xspot... 1x4, x-coord for plotting of model order on time series plots
% yr...	vector of years for the time series, read from col 1 of A 	
% yr1..   vector of valid years for a particular time series
% yspot... 1x4, y-coord for plotting model order on time series plots
% z...	cv, a time series to be analyzed; raw and mean removed

% As set up, the year is the first column of data array A
% How big is this array?

yr = A(:,1);
[nrows,ncols] = size(A);   % how many rows and cols in A
% disp(' nrows  ncols')
% disp([nrows ncols])

%#################  What orders of ARMA models to fit?  #############

% Test only ar(1), ar(2) , ma(1), or arma (1,1)

NN=[1 0
    2 0
    0 1
    1 1];

[m1,n1]=size(ends);
[m,n]=size(NN);
nmods = m;   % how many different order models to test

clc; home;
disp(['A total of ',int2str(nmods),'  models will be tested'])
fprintf('\n');
disp(['AR order    MA order'])
fprintf('\n')
disp(NN)
pause(2)
clc

%################# Pre-allocate some arrays #####################


TBL1=zeros(nmods,5);

parest=zeros(1,8);
ses=zeros(1,8);
mods=zeros(m1,3);
modstat1=zeros(m1,6);
akaik=zeros(nmods,1);
Q=zeros(nmods,1);
pct=zeros(nmods,1);
evar=zeros(nmods,1);
R=zeros(25,nmods);
seax=[0:1:7];
lag=[0:1:24]';

ses(1)=0;
ses(8)=0;
parest(1)=0;
parest(8)=0;
seaxzero=[0 0 0 0 0 0 0 0];


q23df=[35.17 35.17 35.17 35.17]';
q22df=[33.92 33.92 33.92 33.92]';

% ############### A different set of years may be analyzed for each series.

clc, home
opt1 = input('Same rows to be analyzed for all series? Y/N [Y]: ','s');
if isempty(opt1) 
	opt1='Y';
else
end
clc
k5=0;
for k=[colwhich]
k5=k5+1;
%  Optionally use different analysis period for each column.

if opt1=='Y'
	rowgo=1; rowstop=nrows; nyears=nrows;
else
	rowgo=ends(k5,2)-yr(1) +1;
	rowstop=ends(k5,3)-yr(1)+1;
	nyears=rowstop-rowgo+1;
end


% Fill proper year and data into z and yr1

z=A(rowgo:rowstop,k+1);
yr1=yr(rowgo:rowstop);
E=zeros(nyears,5);
E(:,1)=yr1;
% subtract the mean

mz=mean(z);
z = dtrend(z);

%########  COMPUTE ACF OF Z, AND LARGE-LAG 2-STANDARD ERROR ########

ccf=covf(z,25)';     % cross corr fnct, a column vector, lags 0-24
acf=ccf./(ccf(1));     % autocorrelation function

sumr2 = acf(2:25) .* acf(2:25);  % square of acf, lags 1-24 years
cume=cumsum(sumr2);     % cumulative sum of sumr2
cume=cume .*2 + 1;      
sebar= (2/sqrt(nyears))  .* sqrt(cume);

%########## loop here for each model to be tried #################



for k1 = 1:nmods,

nn=NN(k1,:);
nparam = sum(nn)  % total number of ARMA parameters in model


% estimate the parameters of the selected model
% Must call different mfunction depending on whether 0 ar or ma params

if nn(2)==0
	th=ar(z,nn(1));
else
	th = armax(z,nn);
end

% compute the acvf of the residuals of model, and plot acf

[e,r] = resid(z,th); 
E(:,k1+1)=e+mz;
text(0.7,0.8,['Model: (',num2str(nn(1)),',',num2str(nn(2)),')']);
text(0.6,0.9,['Series Number ',num2str(k)]);



% compute the AIC (Salas et al., 1980, p. 97)

% First compute variance of residuals over obs minus the first nparams,

evar(k1)=(std(e(nparam+1:nyears)))^2;   % variance of residuals

akaik(k1)=(nyears-nparam)*log(evar(k1))+2*nparam;



% compute portmanteau statistic based on lags 1-24 of acf of resids

r1= r(2:25)/r(1);
R(2:25,k1)=r1(:);
R(1,k1)=1;
q = length(yr1) * r1*r1';
Q(k1)=q;    % store portmanteau for this model

% store the parameter estimates and their  2 standard errors for later
%  plotting
if k1==1
	parest(2)=th(3,1);
	ses(2)=2*sqrt(th(4,1));
elseif k1==2
	parest(3:4)=th(3,1:2);
	ses(3)=2*sqrt(th(4,1));
	ses(4)=2*sqrt(th(5,2));
elseif k1==3
	parest(5)=th(3,1);
	ses(5)=2*sqrt(th(4,1));
elseif k1==4
	parest(6:7)=th(3,1:2);
	ses(6)=2*sqrt(th(4,1));
	ses(7)=2*sqrt(th(5,2));
end


% compute variance of original series and whitened series, pct var
% of original retained in the whitened series, and fill the table
% describing model results.  Note that r(1) already holds variance of 
% residuals.

varz=(std(z(nparam+1:nyears)))^2
pct(k1)=evar(k1)/varz
pause

TBL1(k1,:)=[nn(1) nn(2) akaik(k1)  q  pct(k1) ];



end   % ##############  end loop of four models for this series #######

k2=0;
modnum=[1 2 3 4]';
while k2~=7
	k2=menu('Select an Item','ACF of Time Series',...
		'50-year Time Plots ','Portmantea & Akaike',...
		'Table','Autocorrs of Resids','Stand Errors of Params',...
		'Ready To Select Model')


	if k2==1    % Plot acf of z
		plot(lag,acf,'+r',lag,zeros(25,1),'w',lag(2:25),sebar,'-b',...
			lag(2:25),-sebar,'-b');
		title(['ACF of Time Series # ',num2str(k)]);
		xlabel('Lag');
		text(.3,.8,'Bars at Two Standard Errors','sc');
		pause
		clg

	elseif k2==2
		% calculate array subscps for plotting at most the most 
		% recent 50 years
			if length(yr1)<=50
				ngo= 1;
			else
				ngo= nyears-49;
			end
		% calculate plotting point for model identifier
			xspot=[0.33 0.83 0.33  0.83];
			yspot=[0.9104  0.9104  0.400  0.400];
		mnp=220;
		for k1=1:4
			mnp=220+k1;
			subplot(mnp)
			plot(yr1(ngo:nyears),z(ngo:nyears)+mz,yr1(ngo:nyears),E(ngo:nyears,k1+1));
			text(xspot(k1),yspot(k1),[num2str(NN(k1,1))  num2str(NN(k1,2))],'sc');
		

			if k1==4
				pause
				clg
			else
			end
		end
	elseif k2==3;  %Plot portmantea and AIC vs model
		mnp=121;
		subplot(mnp)
		plot(modnum,Q,modnum,q23df,'--w',modnum,q22df,':w');
		ylabel('Q Statistic');
		xlabel('Model Number: 1=(1,0)  2=(2,0)  3=(0,1)  4=(1,1)')
		text(1.2,q23df(1),'95% CL, 23 df');
		text(1.2,q22df(1),'95% CL, 22 df');
		mnp=122;
		subplot(mnp);
		plot(modnum,akaik);
		ylabel('AIC');
		pause
		clg
	elseif k2==4
		clc
		disp('        AR        MA       AIC   Porte M.  var(e)/var(z)');
		disp(' ');
		disp(TBL1);
		pause
	elseif k2==5
		clg
		plot(lag,R(:,1),lag,R(:,2),lag,R(:,3),lag,R(:,4),lag,acf,':w');
		xlabel('Lag');
		ylabel('Autocorrelation');
		title('Autocorrelation of Residuals and Original Data');
		text(.5,.8,'red = AR(1)','sc');
		text(.5,.75,'grn = AR(2)','sc');
		text(.5,.70,'blu = MA(1)','sc');
		text(.5,.65,'yel = ARMA(1,1)','sc');
		text(.5,.60,'whi = No Model','sc');
		pause	
		clg
	elseif k2==6
		clg
		plot(seax,parest,seax,seaxzero,'w');
		errorbar(seax,parest,ses);
		title('Parameter Estimates and 2-SE Bars');	
		text(.2,.9,'AR(1)','sc');
		text(.31,.9,'AR(2)','sc');
		text(.42,.9,'AR(2)','sc');
		text(.53,.9,'MA(1)','sc');
		text(.64,.9,'ARMA(1,1)','sc');
		text(.78,.9,'ARMA(1,1)','sc');
		text(.33,.85,'1','sc');
		text(.44,.85,'2','sc');
		text(.66,.85,'AR','sc');
		text(.80,.85,'MA','sc');
		pause
		clg
	elseif k2==7
	end
end

% Select a best model, and store info on it in arrays.

k3=menu('SELECT BEST MODEL','AR(1)','AR(2)','MA(1)','ARMA(1,1)','null');

% Fill in the appropriate info in arrays mods and modstat1 for this series
% These arrays store info on the model fit for poss later use in a table.

if k3<=4 
	mods(k5,:)= [k5 NN(k3,:)];
else
end

modstat1(k5,:)=[0 0 0 0 0 0];

if k3<=4
	modstat1(k5,1:3)=[k5 Q(k3)  pct(k3)]; 
else
	mods(k5,:)=[k5 0 0];
	modstat1(k5,1)=k5;
	modstat1(k5,3)=1.0;
end 

if k3==1
	modstat1(k5,4)=parest(2);
elseif k3==2
	modstat1(k5,4)=parest(3);
	modstat1(k5,5)=parest(4);
elseif k3==3
	modstat1(k5,6)=parest(5);
elseif k3==4
	modstat1(k5,4)=parest(6);
	modstat1(k5,6)=parest(7);
else
end

end

clc
home
disp('You done, Jack!')
disp(' ')
disp('You may want to preserve modstat1 and mods before exiting Matlab')

% pwar.m   prewhitening by AR model, order 1 thru 6


% A... input data array, first col holding year
% e...	residuals of an estimated AR model 
% E...  matrix of residuals for all series
% mz... mean value of original time series z
% ncols... number of columns in A
% nn... lowest AIC model AR order for a particular series
% NN ...  AR orders to be tried (cv).
%	[1 2 3]' says try AR(1), AR(2), AR(3)
% nrows... number of rows in A
%	program assumes A has been trimmed so that the rows contain the
%	valid analysis period.
% opt1... control over period to be analyzed
%	Y- same period (all rows) for all  series
%	N- ends holds beginning and ending years
% q... portmanteau statistic.  See O.D. Anderson, 1976, eqn 9.10, p. 84.
%	Computed on lags 1-24.  Df for test is therefore 24 - AR order.
% TBL1...matrix model results: series #, order model, %var removed, 
%	portmanteau q,
% TBL2...series mean, and estimated parameters and their SE's.
%	Set up now for up to order-3 model only.
% V ... loss functions (top row) for each candidate model
%       model order (row 2)
%       sample size of resids for largest order model (last col)
% varz... variance of z, computed only after deleting number of
%   initial years equal to selected AR order.  This deletion to give
%   better comparison with variance of residuals.
% VMOD... AIC for models tested
% yr...	vector of years for the time series, read from col 1 of A 	
% z...	cv, a time series to be analyzed; raw and mean removed

yr = A(:,1);



[nrows,ncols] = size(A);   % how many rows and cols in A
% disp(' nrows  ncols')
% disp([nrows ncols])



NN=[1 2 3 4 5 ]';
% fprintf('AR orders to be tested:\n')
% disp(NN)

E=zeros(A);
E(:,1)=yr;

TBL1=zeros(ncols-1,6);
TBL2=zeros(ncols-1,8);

% A different set of years may be analyzed for each series.

clc, home
opt1 = input('Same rows to be analyzed for all series? Y/N [Y]: ','s');
if isempty(opt1) 
	opt1='Y';
else
end


% loop thru each of the variables (columns 2,3....ncols-1)

for k=1:ncols


%  Optionally use different analysis period for each column.

if opt1=='Y'
	rowgo=1; rowstop=nrows;
else
	rowgo=ends(k-1,1);
	rowstop=ends(k-1,2);
end


% Fill proper year and data into z and yr1

z=A(rowgo:rowstop,k);
yr1=yr(rowgo:rowstop);

% subtract the mean

mz=mean(z);
z = dtrend(z);

% compute loss function for each order model

V=arxstruc(z,z,NN);
% fprintf('Loss functions for each AR model\n')
% fprintf('Last column is sample size of resids for largest order model\n\n')  
% disp(V);


% select the lowest AIC model

[nn,VMOD] = selstruc(V,'aic');
% fprintf('Lowest AIC AR order = %4.0f\n',nn)
% fprintf('AIC for each model, with model-order below:\n\n')
% disp(VMOD);


% estimate the parameters of the selected model

th = arx(z,nn);


% compute the acvf of the residuals of model, and plot acf

[e,r] = resid(z,th); 
text(0.5,0.9,['Series #',num2str(k-1)])
%pause

% compute portmanteau statistic based on lags 1-24 of acf of resids

r1= r(2:25)/r(1);
q = length(yr1) * r1*r1';


% compute variance of original series and whitened series, pct var
% of original retained in the whitened series, and fill the table
% describing model results.  Note that r(1) already holds variance of 
% residuals.

varz=(std(z(nn+1:length(yr1))))^2;

TBL1(k-1,:)=[k-1  nn  varz  r(1)   (r(1)/varz) q];
TBL2(k-1,:)=[k-1  mz  th(3,1:3)  sqrt(th(4,1:3))];

% replace a col of the zeros matrix E with residuals from model,
% after adding back in the series mean which had previously been 
% removed

E(:,k)=e+mz;

% plot(yr1,z,yr1,e)  
% title(['Time series # ',num2str(k-1),', original-red, whitened-green'])  
% pause
end   % of loop thru each time series


clc
home
disp('    Series     Order    Var(z)    Var(e)       Pct     q-stat' )
fprintf('\n')
disp(TBL1), pause
clc 
home
disp(' Series mean, followed by AR parameters, then their SEs')
fprintf('\n')
disp(TBL2), pause

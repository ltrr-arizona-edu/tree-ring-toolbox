function [I2,stats,c,e,yhat]=stepr(X,y,I1)

%  stepr.m  stepwise regression  

% D. Meko,  10-6-91

% Stepwise multiple linear regression.  Variables entered in order of
% size of correlation with residuals of regression equation containing 
% other variables.

% Variable entry stops when adjusted R squared begins to decrease.

%**************** USER FUNCTIONS NEEDED *****************************

%  sos.m  sum of squares for multiple linear regression


%*****************   INPUT ARGUMENTS   *****************************

% X (m1 x n1)  Matrix from which predictor variables come. [R]
% y (m1 x 1)   predictand variable.  [R]
% I1 (rv)  index to cols of X defining potential predictors. [I]


%****************** RETURNED  **************************************

% I2 (rv)  index to I1 giving variables in final equation.  [I]
% c (cv)  regression coefficientson variables given by I2, constant first. [R]
%		In terms of the initial data array X, the final predictor variables
%		are held in columns  I1(I2)
% e (m1 x 1) residuals from final regression [R]
%		Defined as y-yhat.
% stats (1 x 3)  R-sq, adjusted R-sq, and F-ratio for final equation [R]
% yhat (m1 x 1) predicted values of y from final regression equation [R]		


%*****************  SELECTED OTHER VARIABLES  ***********************************

% c1 - ones cv needed in regression eqn.
% cold - temp vbl holding previous value of c.
% eold - temp vbl holding prev value of e.
% I  - index into R2 of vbl with max correlation.
% I2old - previous step's value of I2.
% I3 - subscripts of X1 giving variables not yet entered.
% k1 - while-loop control. k1=0 until either all vbls are in eqn, or 
%		maximum adjusted R-squared has been reached.
% m1,n1  number of rows,cols in X;  m1 is also number of rows in y.
% m2,n2   number of rows, cols in X1.  X1 is a subset of vbls in X.
% R - highest correlation in R2.
% R1 (n2+1 x n2+1) corr mtx
% R2 (1 X n1) corr coef of e with each col of X1.  Note that e=y before
%		first step in regression.
% statsold - temp vbl holding prev value of stats.
% X1 (m2,n2)  potential predictors 
% X2 (m1,n3) temp mtx of currently selected predictors defined by I2.
% yhatold - temp vbl holding prev value of yhat.


%*********************  NOTES  ********************************************

%  Function assumes you have trimmed X and y to have only valid years, and
%  that the rows of X and y are the years on which the regression is to be
%  done.

%  Variables are added stepwise according to the highest product-moment
%	r between the residuals from the previous step and the variables not
%  yet in the equation.

%**************************************************************************

[m1,n1] = size(X);

X1=X(:,I1);  % Potential predictors
[m2,n2]=size(X1);
e=y;  % Initialize residuals.

I2=[];  %  Start off with no variables  of X1 in the equation.

cold=[0 0];  % Initialize old, or previous, values of coefficients.
statsold=[0 0 0];  % Initialize old values of r-sq, adj r-sq, F-ratio
eold=[];    %  Initialize previous residuals.
yhatold=[];  % Initialize previous predicted values.
I2old=[];  % Initialize previous index to vbls in eqn.

c1=ones(length(y),1);  % Regrn eqn will need a column of ones.

k1=0;    % Loop control for continuing adding variables.
while k1==0,
	R1=corrcoef([e X1]);  %  Correlate residuals with each potential predictor.
	R2 = R1(1,2:n2+1);  % Top row of R1 hold relevant correl coefs.
	[R,I] = max(abs(R2));  % I holds index of largest correl coef in R2 
	I2=[I2 I];  % updated subscripts of predictors from X1
	X2=X1(:,I2);
	[c,stats]=sos(y,X2);
	snew=stats;
	cnew=c;

	if stats(2)<=statsold(2);  % New adj-r-sq less than old.
		stats=statsold;
		c=cold;
		I2=I2old;
		e=eold;
		break;
		yhat=yhatold;
		e=eold;
		k1=1;
	else

		[m3,n3]=size(X2);
		if n3>=n2, k1=1; end  % Exit loop if all vbls have been entered.

		yhat=[c1 X2] * c;
		e=y-yhat;
		
		I2old=I2;   % Store current values of I2,c,stats, yhat, e
		cold=c;
		statsold=stats;
		yhatold=yhat;
		eold=e;

	end;  % Of stats<  loop
end;  % Of while loop

X3=X1(:,I2);





function [I2,I4,stats,c,e,yhat]=stepr2(X,y,I1,delta,num1)
% stepr: stepwise regression
% CALL:  [I2,I4,stats,c,e,yhat]=stepr2(X,y,I1,delta);
%
% D. Meko,  10-6-91;  rev 3-18-92, 1-15-99
%
% Variables entered in order of  size of partial corr between predictand variable 
% and potential  predictors,  removing effect of correlation with vbls already in 
% eqn.   
%
% Criterion for stopping entry of variables depends on number of input args.
%   If 3 input args (no delta), entry stops when adjusted R squared fails to increase
%   If 4 input args, entry stops when R-squared fails to increase by at least delta
%   If 5 input args, entry stops as above according to adjusted R-squared, but with
%     additional constraint that no more than num1 predictors are allowed to enter
%**************** USER FUNCTIONS NEEDED *****************************
%
%  sos.m  sum of squares for multiple linear regression
%  partialr.m 
%
%*****************   INPUT ARGUMENTS   *****************************
%
% X (m1 x n1)  Matrix from which predictor variables come. [R]
% y (m1 x 1)   predictand variable.  [R]
% I1 (rv)  index to cols of X defining potential predictors. [I]
% delta (1 x 1)r  threshold increase in R-squared required for a new
%    predictor to enter equation
% num1 (1 x 1)i  optional.  If there, no more than num1 predictors can enter
%
%
%****************** RETURNED  **************************************
%
% I2 (rv)  index to cols of X giving variables in final equation.  [I]
% I4 (rv) col index saying which of the I1 potential predictors was
%		selected 
% c (cv)  regression coefficientson variables given by I2, constant first. [R]
%		In terms of the initial data array X, the final predictor variables
%		are held in columns  I1(I2)
% e (m1 x 1) residuals from final regression [R]
%		Defined as y-yhat.
% stats (1 x 3)  R-sq, adjusted R-sq, and F-ratio for final equation [R]
% yhat (m1 x 1) predicted values of y from final regression equation [R]		
%
%

%*****************  SELECTED OTHER VARIABLES  ***********************************
%
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

%  First predictor entered has highest pearson-r with predictand. Variables are
%  then added stepwise according to their highest partial-r with the current model
%  residuals (adjusting for removeal of correlation with variables already 
%  predictors in the equation)
%**************************************************************************

[m1,n1] = size(X);

X1=X(:,I1);  % Potential predictors
[m2,n2]=size(X1);
e=y;  % Initialize residuals.

I2=[];  %  Start off with no variables  of X1 in the equation.
I3=1:length(I1);  %  Start off with all variables "not" in equation

cold=[0 0];  % Initialize old, or previous, values of coefficients.
statsold=[0 0 0];  % Initialize old values of r-sq, adj r-sq, F-ratio
eold=y;    %  Initialize previous residuals.
yhatold=[];  % Initialize previous predicted values.
I2old=[];  % Initialize previous index to vbls in eqn.
I3old=I3;  % Initialize pointer to vbls not yet in eqn.


c1=ones(length(y),1);  % Regrn eqn will need a column of ones.

k1=0;    % Loop control for continuing adding variables.  
while k1==0;
	if isempty(I2);  % Special case -- no vbls yet in eqn
	   R1=corrcoef([y X1]);  %  Correlate residuals with each potential predictor.
	   R2 = R1(1,2:n2+1);  % Top row of R1 hold relevant correl coefs.
      [R,I] = max(abs(R2));  % I holds index of largest correl coef in R2 
      I2add=I3(I); % will later add this variable to set now in equation
	   I3(I)=[];
	else;  % At least 1 vbl has been previously entered
      R1=partialr([y X1(:,I3)],X1(:,I2));
      R2=R1(2:length(R1));
      [R,I] = max(abs(R2));
      I2add=I3(I);
	   I3(I)=[];  % remove the newly used predictor from list of not in
  end
	I2=[I2 I2add];  % updated subscripts of predictors from X1 now included in eqn
	X2=X1(:,I2);
	[c,stats]=sos(y,X2);
	snew=stats;
	cnew=c;
   
   if nargin==3;
      criter=stats(2);
      criterold = statsold(2);
      npred=n2+1; % make this large enough to cover all potential predictors
   elseif nargin==4;
      criter = stats(1)-delta;
      criterold = statsold(1);
      npred=n2+1; % see above
   elseif nargin==5;
      criter=stats(2);
      criterold=statsold(2);
      npred=num1; % you have specified num1 as max allowable number of predictors
   else
      error('number of input args must be 3,4 or 5');
   end
   
   if ((criter<=criterold) & (length(I2)>1)) | (length(I2)>1 & length(I2old)==npred);
      % Criterion on R square or adj-r-sq satisfied and the variable just entered 
      % is not the first to enter -- or -- variable just entered is not the first
      % to enter and number of predictors already in model at last step is npred.
      k1=1; % will not want to go thru while k1 again
      % Set values of stats,c,I2,e,yhat to those for previously fit model.
      % In other words, override the most recent results because no gain from model
      % expansion.
		stats=statsold;
		c=cold;
		I2=I2old;
		I3=I3old;
		e=eold;
		yhat=yhatold;
	else; %OK, some gain in adj R-sqrd from the new model. So get yhat and e
      
		yhat=[c1 X2] * c;
      e=y-yhat;
      
      % Will want to exit the while loop if all variables have been entered
      [m3,n3]=size(X2);
      if n3>=n2, 
         k1=1; 
         break;
      else; % will keep going on, so store some values for current model
         I2old=I2;   % Store current values of I2,c,stats, yhat, e
         I3old=I3 ;  % Store...
         cold=c;
         statsold=stats;
         yhatold=yhat;
         eold=e;
      end; % of if n3>=n2
      
	end;  % Of stats<  loop
end;  % Of while loop

I4=I2; % relative predictor  col index specific to the set of variables X(I1)
I2=I1(I2); % corresponding col index pointing to original matrix X

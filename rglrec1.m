function [c,IW,S,yhat]=rglrec1(X,Y,IY,IX,nX,yrs)

% Gridpoint reconstruction of pdsi.  Spatial weighting.

% D. Meko  11-23-92


%****************  USER-WRITTEN FUNCTIONS CALLED DIRECTLY
% 
% stepr1.m  stepwise regression
% 


%****************  INPUT ARGS
%
% X (m1 x n1)   potential predictor data, each row a year, col a variable
% Y (m2 x n2)   predictand data (potential)
% IY (m3 x 1)   index to cols of Y telling which vbls to use for pdctands
% IX (m5 x n5)  index to cols of X telling which subset of variables are
%		set of potential predictors.  m5=m3.
% nX (m4 x 1) number of predictor variables in each model corresp to rows
%		of IY and IX.  m4=m3.
% yrs(4x2)  beginning and ending years of :
%	row1  X
%	row2  Y
%	row3  calibration period
%	row4  long-term desired reconstruction


%*****************  OUTPUT ARGS
%
% c (m3 x (n5+1))  estimated regr coefficients, constant first, remainder
%		in order as specified by entries in IX.  Zero-filled if model
%		smaller than n5 vbls.
% IW (m3 x n5)  index to cols of X specifying which variables were
%		were selected as predictors.  Zero-filled as needed.
% S (m3x3)  R-squared, adjusted R-squared, F-ration for model
% yhat (m6 x m3) long-term reconstruction, rows corresp to years
%		in yrs(4,:)


%******************  NOTES
%
% Variables entered stepwise according to max correlation with residuals
% from current model.  Details in doc for stepr1.m.
%
% First use:  I had generated 155 gridpoint values of pdsi, and 250
%  single-site reconstructions of pdsi at tree-ring sites.  Needed to
%  weight the single-site pdsi reconstructions into gridpoint 
%  reconstructions.  Also had run functions to get appropriate indices
%  of tree-ring sites in search radius for each gridpoint.

%******************  YEAR POINTERS

yr1= (yrs(1,1):yrs(1,2))';  % cv of years of X
yr2= (yrs(2,1):yrs(2,2))';  % ... of Y
yr3= (yrs(3,1):yrs(3,2))';  % ... of calib period
yr4= (yrs(4,1):yrs(4,2))';  % ... of long-term recon

LX1= yr1>=yrs(3,1) & yr1<=yrs(3,2);  % pointer to calib-per years in X
LX2= yr1>=yrs(4,1) & yr1<=yrs(4,2);  % pointer to long-term rec in in X
LY1= yr2>=yrs(3,1) & yr2<=yrs(3,2);  % pointer to calib-per years in Y


%******************   PREALLOCATE, SIZE, INITIALIZE

[m1,n1]=size(X);
[m2,n2]=size(Y);
m3=length(IY);
m4=length(nX);
[m5,n5]=size(IX);
yhh=length(yr3);

a=NaN;
m6=length(yr4);
yhat=a(ones(m6,1),ones(m3,1));  
S=a(ones(m3,1),ones(3,1));
c=zeros(m3,n5+1);
IW=zeros(m3,n5);
cones= ones(length(yr4),1);   % col vect of ones for regression 



%****************   REGRESS AND RECONSTRUCT

for i=1:m3;  % Loop for each predictand
	i
	y=Y(LY1,IY(i));  % pull the correct rows and col for predictand
	I1= IX(i,1:nX(i));  % a rv specifying potential predictors

	[I2,stats,cc,ee,yh]=stepr2(X(LX1,:),y,I1);
	IW(i,1:length(I2)) = I2;
	yhat(:,i)=[cones X(LX2,I2)] * cc;  % generate long-term reconstruction
	S(i,:)=stats;
	c(i,1:length(cc))=cc';

	%plot(yr4,yhat,yr3,yh);  % plot reconstructed, overlaying for a 
		% double check the predicted values for calib period.

	%pause
end

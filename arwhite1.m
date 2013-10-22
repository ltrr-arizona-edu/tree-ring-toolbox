function [K1,VR,COEF,ER,E]=arwhiten(xdata)
% arwhiten.m   .... sample code to illustrate use of function
%     whit1.m .  Prewhiten cols of a matrix and gives
%     a graphic display of: 

%  1. acf of residuals from ar model, with 99% confidence bands
%  2. time series plots of original and whitened data, full length
%  3. Last nshort years of the same time series 

% By D. Meko 3-23-92

%******  USER-WRITTEN FUNCTIONS NEEDED
%
%  whit1.m 

%*********  MATRIX OUTPUT
%
%	E (m,n) AR residuals of each series (same size as input data matrix
%	VR (n x 1)  ratio of variance of residuals to orig series
%	K1 (n x 1) selected order of AR model
%   COEF (n x 6)  AR coefficients each series ;  non-applicable cols
%		zero filled
%   ER (n x 6) two-standard error limits around coefficients in COEF
%   


% Read the descriptions of arguments in the function comments to 
% acf.m, pacf.m , whit1.m  before using this code.


k2=input('DOES COL 1 CONTAIN THE YEAR? Y/N [N] ','s');
if isempty(k2), k2='N';  end;
if k2=='Y',
	xyears=xdata(:,1);
	xdata(:,1)=[]; % get rid of the year column
else
	yr=input('START, END YEARS FOR DATA MATRIX:  ');
	xyears=(yr(1):yr(2))';
end;
nshort=input('HOW MANY YEARS FOR CLOSE-UP VIEW OF RESIDUALS? ');


[m,n]=size(xdata);   % Compute number of rows, cols in xdata*

COEF=zeros(n,6);  % Preallocate, allowing room for 6 coefs
ER=zeros(n,6);  % Will hold 2-std errors for coefs
VR=zeros(n,1);  %  Will hold ratios of resid to orig variance
K1=zeros(n,1);  % Will hold selected order of ar model
E=zeros(m,n);   %  Will hold ar residuals (whitened series)

for i=1:n;  % Loop for each col of xdata (each time series)
	[e,k1,vrat,arcs]=whit1(xdata(:,i),5,1);
	
	text(4,.2,['SERIES # ',int2str(i)]);
	text(1,.1,['AR ORDER = ',int2str(k1)]);
	
	pause

	E(:,i)=e;   % store residuals of each series as cols of E
	VR(i)=vrat;  % store ratios of variance of residuals  to orig series
	K1(i)=k1;  % store selected order of AR model
	
	if k1 ~=0 ;  % if null model has not been selected, store coefs and
		% their 2-standard-errors
		COEF(i,1:k1)=arcs(1,1:k1);
		ER(i,1:k1)=arcs(2,1:k1);
	end

	clg
	subplot(211)
	plot(xyears,xdata(:,i),'-',xyears,E(:,i),'--');
	title(['ORIG-SOLID, WHITENED-DASHED, SERIES ',int2str(i)])
	pause
	y1=length(xyears)-(nshort-1);
	y2=length(xyears);
	yr1=(xyears(y1):xyears(y2));
	plot(yr1,xdata(y1:y2,i),'-',yr1,E(y1:y2,i),'--');
	pause
	subplot
 end
	
if k2=='Y';  % put a year col back in residual array E, if data array 
%              had a year col
	E=[xyears E];
end

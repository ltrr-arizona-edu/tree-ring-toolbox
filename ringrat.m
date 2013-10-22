function Y=ringrat(X)
% Ratio of ring-width change to previous rw, for use in cross dating


%*******  INPUT ARGS
%
% X (mX x nX) nX time series of ring width, all same length


%******* OUTPUT ARGS
%
%  Y (mY x nY) corrsponding difference series
%	mY will be one less than mX



%***  NOTES
%
% To avoid infinite values, if t-1 is missing, function sets denominator of ratio to 
% one-half the smallest non-zero measurment in the data.
%
% If both ring t and ring t-1 are missing, ratio is set equal to 1.


[mX,nX]=size(X);


% Form the time t and time t-1 arrays
X1 = X(2:mX,:);  % time t
X2 = X(1:mX-1,:); % time t-1

% Form logical pointer to denominator values not equal to zero
% Calc one half the minimum non-zero value for each series
LLa=X1>0;
for i=1:nX;
	colx=X(:,i);
	LL=LLa(:,i);
	xx(i)=min(colx(LL))* 0.5;

end
Y=xx;

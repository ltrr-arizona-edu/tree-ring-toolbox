function [fsig,fuse]=ftest(f95,df1,df2);

% Given numerator degrees of freedom (df1) and denominator
% degrees of freedom (df2), returns the points of the F
% distribution for 95% significance (fsig) and for practical 
% usefulness as predictor (fuse).  

% Requires that mat file ftable95.mat be loaded beforehand.

% Handles df1,df2 either scalar or vector (rv)

% fuse comes from Draper and Smith (1981, p. 93), based on 
% J. M. Wetz 1964 phd thesis.  Wetz finds that for the range of
% predicted values to be practically large compared to the standard
% error of the response variable, need an F-ratio at least 4 times
% the selected percentage point of the F-distribution. 

% So say testing a correlation coef.  Find r=0.5 for 250 years of data.
% Appropriate mean square ratio (Panofsky and Brier 1968, p. 93) is
%   [(rsquared / 1)]  /  [(1-rsquared)/(N-2)], or
%   [.25 / 1] /  [.75 /248]  = 248/3 = 82.67

% Call to this routine, with appropriate values
%   df1=1
%   df2=248
% gives  fsig=3.88 and fuse=15.53.  

% The F-ratio is therefore both statistically signif and practically 
% useful at the 95% level.


nsize=length(df1);  % can do a vector of tests in one call
fsig=zeros(1,nsize);

if nsize == 1
	fsig=table2(f95,df2,df1);
else
	for i=1:nsize
		fsig(i)=table2(f95,df2(i),df1(i));
	end
end

fuse=4*fsig;

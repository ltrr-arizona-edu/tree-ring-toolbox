function [y,iy,b]= filtbin(x,ix,m,kopt)
% filtbin: binomial filtering 
% CALL:  [y,iy,b]= filtbin(x,ix,m,kopt);
%
% Meko 4-4-98
%
%************** IN
%
% x (mx x 1)r  time series, mx observations
% ix (mx x 1)r  time vector corresponding to x (e.g., (1:366)' for daily means)
% m (1 x 1)i  desired length of binomial filter; must be odd and 3 or larger
% kopt (1 x 1)i  options
%   kopt(1)  mode for handling ends of time series
%      ==1 no augmentation of series before filtering
%          (output series is shorter than input because of loss of end years)
%      ==2 circular: assumes input series repeats itself over and over 
%      ==3 extend ends with long-term means before filtering
%
%********** OUT 
%
% y (my x 1)r filtered time series, my observations
% iy (my x 1)r time vector for y
% b (1 x mb)r binomial filter weights
%
%************ NOTES 
%
% First application was to smooth long-term daily mean geopotential height data. Had
% 366 values corresponding to days of year.  For this application, circular option 
% kopt(1)==2 is best.
%
% Filtered value at time iy(i) is centered value resulting from weighting
% earlier and later values
%
% Number of filter weights must be odd;  even number will make function bomb

%---------  SIZE
[mx,nx]=size(x);
[mtemp,ntemp]=size(ix);
if mx~=mtemp | nx~=1 | ntemp~=1;
   error('x and ix must be col vectors of same length');
end
nobs1 = mx;
clear mtemp ntemp mx nx 

[mtemp,ntemp]=size(m);
if mtemp~=1  | ntemp~=1;
   error('m must be scalar');
end
if m>nobs1 | m<3;
   error('m must be greater than 2 and smaller than length of series c');
end
if mod(m,2)==0;
   error('m must be odd');
end


%***************  COMPUTE FILTER WEIGHTS

b=[0.5 0.5];
kwh1 = 1;
while kwh1==1;
   b=conv([0.5 0.5],b);
   mb = length(b);
   if mb == m;
      kwh1=0;
   else
   end
end
nshift=fix(m/2);

%****************  EXTEND ENDS OF TIME SERIES
xmn = mean(x);
iy = ix;
if kopt(1)==1;
   xaug = repmat(NaN,nobs1,1);
   iy((nobs1-nshift+1):nobs1)=[];
   iy(1:nshift)=[];
   
elseif kopt(1)==2;
   xaug = x;
elseif kopt(1)==3;
   xaug =repmat(xmn,nobs1,1);
else
   error('kopt(1) must be 1,2, or 3');
end

z = [xaug; x; xaug];


%***************** FILTER TIME SERIES
y = filter(b,1,z);
y = y((nobs1+nshift+1):(2*nobs1+nshift));
Lnan = isnan(y);
y(Lnan)=[];




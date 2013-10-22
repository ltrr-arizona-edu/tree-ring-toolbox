function Z=clicon2(X)
% Z=clicon2(X)
% 
% Convert EarthInfo style monthly ppt data  
% into format used by Meko functions
%
%********************* IN and OUT ARGS *************************************
%
% X (mX x nX)i monthly ppt data keypunched from GS64.  In inches, 
% with -9.99 as missing value code
%
% Y (mY x nY)r same data, but  with
%  NaN as missing value code

a=NaN;
close all

[mX,nX]=size(X);
yr=X(:,1);
yrgo=yr(1);
yrsp=yr(mX);
if (yrsp-yrgo+1) ~= mX
	error('Row size of X inconsistent with start, end years');
end

% Display start, end years
str1=sprintf('Years: %4.0f-%4.0f\n',yrgo,yrsp)

% Relace -9.9900 with NaN
w=X(:,2:13);
L=w==-9.9900;
s=sum(sum(L));
if s>0;
	w(L)=a(ones(s,1),:);
end


% Slap on year
Z=[yr w];
 
% Boxplot of monthly data
figure(1);
boxplot(w)
xlabel('Month')
ylabel('PPT (inches)')

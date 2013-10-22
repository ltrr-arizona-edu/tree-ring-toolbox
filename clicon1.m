function Z=clicon1(X)
% Z=clicon1(X)
% 
% Convert spreadsheet-made monthly ppt data from Green and Sellers (1964) 
% into format used by Meko functions
%
%********************* IN and OUT ARGS *************************************
%
% X (mX x nX)i monthly ppt data keypunched from GS64.  In hundredths of
%	inches, with -9 as missing value code
%
% Y (mY x nY)r same data, but in inches (with decimal point), and with
%  NaN as missing value code

a=NaN;

[mX,nX]=size(X);
yr=X(:,1);
yrgo=yr(1);
yrsp=yr(mX);
if (yrsp-yrgo+1) ~= mX
	error('Row size of X inconsistent with start, end years');
end

% Display start, end years
str1=sprintf('Years: %4.0f-%4.0f\n',yrgo,yrsp)

% Relace -9 with NaN
w=X(:,2:13);
L=w==-9;
s=sum(sum(L));
if s>0;
	w(L)=a(ones(s,1),:);
end

% Divide by 100 to convert to inches
w=w/100;

% Slap on year
Z=[yr w];
 
% Boxplot of monthly data
figure(1);
boxplot(w)
xlabel('Month')
ylabel('PPT (inches)')

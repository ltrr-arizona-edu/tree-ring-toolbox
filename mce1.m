function [S,hmax]=mce1(X,yrs1,yrs2)
% mce1:  approximation to minimum covering ellipsoid for identifying extrapolations
% CALL: [S,hmax]=mce1(X,yrs1,yrs2);
%
%************  INPUT ARGS **********************
%
% X (m x n) predictor matrix, ones included as col 1
%      m years, n columns; includes calib data and long-term data
% yrs1 (1 x 2) years covered by X
% yrs2 (1 x 2) years in calibration period
%
%************  OUTPUT ARGS  ******************
%
% S (m x 3)  statistical summary showing extraps and interpols
%    col 1 - year
%    col 2 - h measure (Weisberg's h-star on p. 229 and 237
%    col 3 - 1 for extrap, 0 for interpolation
% hmax  (1 x 1) maximum value of hstar for calibration period
%	Extrapolation is defined as hstar>hmax
%
%************* NOTES 
%
% Source- Weisberg (1985), p. 237.  Note that this
% is not the smallest -volume ellipsoid of Titterington (1978)

[m,n]=size(X);


yr1=[yrs1(1):yrs1(2)]';  % cvs of years in full and calib
yr2=[yrs2(1):yrs2(2)]';

hc=zeros(length(yr2),1);
h=zeros(m,1);
S=zeros (m,3);

L=yr1 >= yrs2(1) & yr1 <= yrs2(2);  % pointer to calib-pd rows

Z=X(L,:); % pull calib sub-data out of X.

II = inv(Z' * Z);  % as in Weisberg, but Z for X

for i = 1:length(yr2)
	hc(i)= Z(i,:) * II * (Z(i,:))';  % construction sample h
end

hmax=max(hc);

for i=1:m
	h(i)= X(i,:) * II * (X(i,:))';
end


LL=h>hmax;

S=[yr1,h,LL];

%plot(S(:,1),S(:,2),S(:,1),hmax(ones(m,1),:))
%ylabel('h-star');  xlabel('Year');
%title('h-star distances, from Weisberg 1985, 237.  From mce1.m')
%pltext(.1,.9,'Values above line are extrapolations')
%pause


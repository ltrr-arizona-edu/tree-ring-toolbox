function Y=dec2dms(X)
% lat-long file in decimal degrees to degrees, min seconds


%*******   INPUT ARGS
%
% X (mX,2) - decimal degrees of lat, long for mX locations 

[mX,nX]=size(X);

XD=fix(X);  % round toward zero for degrees
XM1=(X-XD)*60;

k1=input('Round to nearest minute? Y/N [Y]','s')
if(isempty(k1), k1='Y'; end;
if(k1=='Y' | k1=='y')
	XM=round(XM1);  % to nearest whole degree
	XS=zeros(mX,2);  % seconds are zero
else
	XM=fix(XM1);
	XS1=(XM1-XM)*60;
	XS=fix(XS1);
end
 


Y=[XD(:,1) XM(:,1) XS(:,1)  XD(:,2) XM(:,2) XS(:,2)];


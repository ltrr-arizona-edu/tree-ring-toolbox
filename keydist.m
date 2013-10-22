function d=keydist(A,m);
% keydist: subfunction of regcli3a.m giving threshold distance for mth 
%  nearest station in matrix of distances of points to stations
%
%***********8  IN ********
%
% A (mA x nA)r distances from mA points to nA stations
% m (1 x 1)i  desire distance to the mth nearest station for each point
%
%*********** OUT
%
% d (mA x 1)r  distance to the mth nearest station for each point


[mA,nA]=size(A);
B=A';
if any(any(isnan(B)));
   error('Function as written does not allow NaN in A');
end

if m>nA;
   error('m must be equal to or smaller than col dim of A');
end


d=repmat(NaN,mA,1);


[Y,I]=sort(B);
i = I(m,:);
j = sub2ind([nA mA],i,[1:mA]);
b=B(:);
d = b(j);

   



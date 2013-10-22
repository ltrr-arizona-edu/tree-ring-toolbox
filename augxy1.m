function B=augxy1(A,i,what,akey)
% augxy1:  augment lat & long vector by 1 by interpolation
% CALL:  B=augxy1(A,i,what,akey);
%
% Meko 9-6-97
%
%******************  IN 
%
% A (? x 2)r  matrix of latitude (col 1) and longitude (col 2)
% i (mi x 1)i  rows of A defining segment over which interp to happen
% what (1 x ?)s  which variable to interpolate from (long or lat)
% akey (1 x 1)r  the value of that variable at whic other variable at which estimate needed
%
%*****************  OUT 
%
% B (? x 2)r augmented matrix of lats and longs
%
%
%*************** NOTES 
%
% i must be monotonic increasing
% what must be 'lat' or 'lon'
%
%
% Example:
% B=augxy(A,[4 5 6 7]','long',34)


if ~all(diff(i)==1);
   error('in arg i must be monotonic increasing by 1');
end


[mA,nA]=size(A);

% pull segment to work on
A1 = A(i,:);
nin = length(i);
lat0=A1(:,1);
lon0=A1(:,2);

Atop=A(1:(i(1)-1),:);
Abot=A((i(nin)+1):mA,:);


if strcmp(what,'lat'); %  LONGITUDE NEEDS TO BE ESTIMATED, GIVEN LAT
   xlon = intrplon(lat0,lon0,akey);
   
   % find which row to slip new value into
   L1 = lat0<akey;
   L2 = lat0>akey;
   
   if  all(diff(lat0)>=0);  % A1 lats were increasing down the row
      B1= [A1(L1,:); [akey xlon]; A1(L2,:)];
   elseif all(diff(lat0)<=0);  % A1 lats were decreasing down the row
      B1= [A1(L2,:);  [akey xlon];  A1(L1,:)];
   else
      error('Input lon segment must be monotonic');
   end
elseif strcmp(what,'lon'); % LAT NEEDS EST, GIVEN LONG
   xlat = intrplat(lon0,lat0,akey);
   
   % find which row to slip new value into
   L1 = lon0<akey;
   L2 = lon0>akey;
   
   if  all(diff(lon0)>=0);  % A1 longs were increasing down the row
      B1= [A1(L1,:); [xlat,akey]; A1(L2,:)];
   elseif all(diff(lon0)<=0);  % A1 longs were decreasing down the row
      B1= [A1(L2,:);  [xlat,akey];  A1(L1,:)];
   else
      error('Input lon segment must be monotonic');
   end
end

B=[Atop; B1; Abot];   
   
   
      




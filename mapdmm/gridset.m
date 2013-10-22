function X=gridset(lat1,lon1,inc);
% gridset:  build a lat/lon grid whose edges are defined by lats & lons of points
% X=gridset(lat1,lon1,inc);
% Last revised 2-16-01
%
% Needed to build grid for Paclim 2001 mtg.  Have locations of 60 pcp stns in SE Arizona.  Want to 
% use grid system to define "wetness" of summers using daily P data
%
%*** INPUT
%
% lat1 (? x 1)r  latitudes of points (dec deg)
% lon1 (? x 1)r  longitudes of points, W long negative
% inc (1 x 1)r   desired increment of lat and lon :  spacing between gridpoints
%
%*** OUTPUT
%
% X(? x 2)r  lat, lon of gridpoints, arranged S to N, W to E
%
%*** REFERENCES
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE


% Get min and max of point lat & lon
latmin = min(lat1);
latmax=max(lat1);

lonmin=min(lon1);
lonmax=max(lon1);

%  Initialize grid to extend beyond points
lata=[floor(latmin):inc:ceil(latmax)]';
lona=[floor(lonmin):inc:ceil(lonmax)]';


% Shrink grid inside of point extremes


L1 = lata>=latmin & lata<=latmax;
lat2=lata(L1);

L1 = lona>=lonmin & lona<=lonmax;
lon2=lona(L1);

nlat=length(lat2);
nlon=length(lon2);

a = lat2'; % rv of latitudes
A= repmat(a,nlon,1); % dupe rows
a = A(:);

b = lon2;
B = repmat(b,1,nlat);
b = B(:);

X = [a b];

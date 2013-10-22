function [latrng, lonrng]=xyedge(A,distkm)
% xyedge: latitude and longitude
% [lat lon]=xyedge(A,distkm); boundary of region enclosing specified points
% Last revised 2-10-01
%
% You have a network of tree-ring sites and want to find all climate stations within a specified 
% search distance of any of the sites.  Also to include any sites in the lat/long box defined by
% the outer latitudes and longitudes of sites.  To do this, you need a lat/lon rectangle to search your 
% climate database with.  xyedge.m computes the lat and lon limits of the required box
%
%
%*** INPUT
%
% A (mA x 2)r  latitude and longitude (dec deg, west negative long) of tree-ring sites
% distkm (1 x 1)r search radius (km) around each site (see notes)
%
%
%*** OUTPUT
%
% latrng (1 x 2)r lower and upper latitudes of box enclosing sites + search radius
% lonrng (1 x 2)r lower and upper longitudes of box
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED
%
% gcdist
%
%*** TOOLBOXES NEEDED -- MAPPING
%
%*** NOTES
%
% latrng, lonrng define a box that enclose the sites plus the specified search radius around each

% 

% Convert search distance to arc degrees
distdeg = km2deg(distkm);


[mA,nA]=size(A); % mA is number of sites in set


% Allocate for merged lat ,lon vectors
Blat=repmat(NaN,100*mA,1);
Blon=repmat(NaN,100*mA,1);

for n = 1:mA; % loop over points
    i1=1+(n-1)*100;
    i2=i1+99;
    a = A(n,:);
    lat=a(1);
    lon=a(2);
    [latc,lonc] = scircle1(lat,lon,distdeg); % gives 100 points
    Blat(i1:i2)=latc;
    Blon(i1:i2)=lonc;
end;

latrng=[min(Blat) max(Blat)];
lonrng=[min(Blon) max(Blon)];
function gage (lat,lon,mplatlim,mplonlim,scale,color)
% gage: add symbol for a stream gage to a map
% CALL: gage (lat,lon,mplatlim,mplonlim,scale,color);
%
% Meko 8-18-98
%
%****************  IN 
%
% lat(1 x 1)r  latitude of center of gage 
% lon(1 x 1)r  longitude ...
% mplatlim (1 x 2)r  MapLonLimit propert of map
% mplonlim (1 x 2)r  MapLonLimit property of map
% scale (1 x 1)r  decimal fraction of lat range of map for width of gage symbol
%     (e.g., 0.02)
% color(1 x 3)r  fill color for the gage, which is a patch object (e.g., [1 0 0]
%
%************** OUT
%
% No output args.
% Action is to add an "upside down triangle" patch object to a map
%
%************** NOTES
%
% User might need to vary scale till satisfied with relative size of gage 
% symbol.  I have been happy with scale==0.02 on one map.



% Calculate width and heigth of the symbol in degrees
delx = diff(mplonlim) * scale;
dely = (diff(mplatlim)/diff(mplonlim)) * delx;

% Compute x, y plotting points for corners of symbol
% Points 1, 2, 3 correspond to upper left, upper right, and bottom center of gage 
x1=lon - delx/2;
x2=lon + delx/2;
x3=lon;

y1 = lat + dely/2;
y2= lat + dely/2;
y3=lat - dely/2;

% Set the latlon for the symbol
xy = [y1 x1; y2 x2; y3 x3];

% plot the symbol
patchm(xy(:,1),xy(:,2),color);

   

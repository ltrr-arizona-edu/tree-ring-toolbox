function x=scalebar(neglon,lat,dist)
% x=scalebar(neglon,lat,dist)
% Get rightmost (eastern edge) negative longitude value for use in 
% plotting a scale bar on a surfer map
%
% D Meko 12-6-96
%
%******************* IN ARGS ************************************
%
% neglon (1 x 1)r  negative longitude (e.g., -124) for western
%		edge of bar, decimal degrees
% lat (1 x 1)r representative latitude, in decimal degrees, that
%		scale will apply to
% dist (1 x 1)r   distance, in km, that you want the scale bar to
%	cover (e.g., 100)
%
%
%***************** OUT ARG ****************************************
%
% x (1 x 1)r negative longitude  for right side of bar 
%		(e.g., -122.8)
%
%
%************* NOTE *****************************
%
% This simpleton algorithm starts with a bar of zero length in longit units,
% then repeatedly increments the longitude of west edge by 0.01 until
% distance from left to right edge exceeds dist.  Thus, might go ...
%  -124, -123.99, -123.98,....
%
% This accuracy should be ok for most maps.  Only a problem if wholw
% map covers, say, only a tenth of a degree



p1=[neglon lat];
p2init=[neglon lat];
x=neglon;
d=0;


p2=p2init;
while (dist-d)>0.01*dist;
	p2=[x lat];
	d=gcdist(p1,p2);
	x=x+0.01;
end

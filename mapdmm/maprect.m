function maprect(latlim,lonlim)
% maprect:  draw a rectangle on a map to allow identifying where to extend lines to an internal map
% maprect(latlim,lonlim);
% Last revised:  3-6-01
%
% State, country, etc;, lines may not extend to lat and lon limits; need to extend them
%
%*** INPUT
%
% latlim (1 x 2)r  latitude limits of the new, inner map
% lonlim (1 x 2)r  longitude limits...
%
%
%*** OUTPUT
%
% No arguments
% The rectangle defined by latlim, lonlim is drawn on the on-screen map.  User can then see where
% existing state lines intersect the inner boundary.  Can zoom and use inputm to record those points. 
% Then can modify map script to extend lines if needed to the boundary
%
%*** REFERENCES
%*** UW FUNCTIONS CALLED -- NONE
%
%*** TOOLBOXES NEEDED -- MAPPING
%
%*** NOTES
%
% Assumed that a map is in the current figure window when rin maprect


lat = [latlim(1)  latlim(1)  latlim(2)  latlim(2) latlim(1)];
lon = [lonlim(1)  lonlim(2)  lonlim(2)  lonlim(1) lonlim(1)];

htestline = linem(lat,lon);
set(htestline,'Color',[1 0 0]);
function [lat1,lon1]=maptrim1(lat0,lon0,latlim,lonlim,maxdiff)
% maptrim1:  adjusts long and lat data so that lines come up to trimboundary
% CALL: [lato1,lon1]=maptrim1(lat0,lon0,latlim,lonlim,maxdiff);
%
%**************** IN  
%
% lat0 (m1 x 1)r  latitude; separated by NaNs if segments
% lon0 (m1 x 1)r  longitude corresp to lato
% latlim (1 x 2)r  latitudes of boundaries of trim region
% lonlim (1 x 2)r  longitudes of boundaries of trim region
% maxdiff(1 x 1)r maximum desired separation (angular degrees) between any two points on 
%   output lat1,lon1
%
%************** OUT 
%
% lat1 (m2 x 1)r  augmented (possibly) latitudes
% lon1 (m2 x 1)r  longitude corresponding to lat1
%
%
%*********** NOTE 
%
% Need:  use of MATLAB function maptriml sometimes results in huge gap in line
% when line should go to trim boundary.  This is most obvious if the line is a 
% straight line segment (e.g., the S Colorado border).  In that case, maybe only
% two x,y pairs define the line.  Trimming results in dropping off the whole line.
%
% Approach:  Loop over the line segments.  Look for segments that cross a boundary.
% For such a segment, linearly interpolate so that all points on the line are no more
% than maxdiff apart
%
%

%**************** SIZE

[m1,n1]=size(lat0); % m1 is number of points
[mtemp,ntemp]=size(lon0);
if n1~=1 | m1 ~= mtemp;
   error('lat0 and lon0 must be col vectors of same length');
end

[m3,n3]=size(latlim);
[m4,n4]=size(lonlim);
if m3~=1 | m4 ~=1 | n3~=2 | n4~=2;
   error('latlim and lonlim must be 1 x 2');
end
clear m3 n3 m4 n4 mtemp ntemp
lon1=[];
lat1=[];

%****************** COMPUTE NUMBER OF LINE SEGMENTS
%
% Seems that the mapping toolboxes lat and long data files sometimes start with a NaN
% and sometimes end with a NaN.  The number of NaNs should equal the number of line
% segments.  I try to make things consistent by making all files start with NaN
% rather than end with a NaN

% Count line segments
Lnan = isnan(lat0);

% Adjust endpoint NaN location
if isnan(lat0(1)) & isnan(lat0(m1));
   lat0(m1)=[];
   lon0(m1)=[];
   m1=m1-1;
end
if ~isnan(lat0(1)) & ~isnan(lat0(m1));
   error('No NaN at either end of lat0');
end
if isnan(lat0(1));
   % No action needed; NaN is first value
elseif isnan(lat0(m1));
   lat0(m1)=[];
   lon0(m1)=[];
   lat0=[NaN;  lat0];
   lon0=[NaN; lon0];
else
   error('No comprendo this one');
end

%**************** COMPUTE START AND END ROW FOR EACH SEGMENT
%
% sements begin with NaN; have nseg segments

Lnan = isnan(lat0);  % revise the locations of NaNs
% check that longs have NaNs is same spots as lats
nseg = sum(Lnan);

Ltemp = isnan(lon0);
Ltemp = Ltemp & Lnan;
if sum(Ltemp)~=nseg | sum(Ltemp)~=nseg;
   error('NaNs in longs and lats not in identical rows');
end

i1 = find(Lnan);  % rows with NaNs
% Compute starting rows of each segment
igo = i1+1;

% Compute ending rows
i2 = i1;
i2(1)=[]; 
isp = [(i2-1);  m1];

% Compute number of points in each segment
npoints = isp - igo +1;
if min(npoints)<2;
   error('A line segment has fewer than 2 defining points');
end


%************************ AUGMENT THE LAT AND LONG FILES
lat1=[];
lon1=[];

% Loop over the segments
for n = 1:nseg;
   i1go = igo(n);
   i1sp = isp(n);
   
   % Get the segments
   x = lon0(i1go:i1sp);
   y = lat0(i1go:i1sp);

   
   % Check whether segment touches inside the lat and long limits
   Lin  =  x>lonlim(1) & x<lonlim(2) & y>latlim(1) & y<latlim(2);
   if all(Lin) | ~any(Lin); % no need to adjust; points all inside or outside bry
      lat1=[lat1; NaN; y];
      lon1=[lon1; NaN; x];
      
   else; % line touches or crosses boundary
      [y1,x1]=interpm(y,x,maxdiff);
      nx = length(x1)-length(x);
      lat1 = [lat1; NaN; y1];
      lon1 = [lon1; NaN; x1];
   end
end


disp('here');





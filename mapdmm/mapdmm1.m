function mapdmm1(coast)
% mapdmm1: make a map out of a usgs-formatted matlab coast file
%CALL: mapdmm1(coast,w,s,c);
%
% Meko 7-14-97
%
%************************************************************
%
% coast (? x 2)r  long, lat for points along boundary, signed long/lat units
%   Line segments begin with a line of two NaN's
% s,w,c ( all string) line style, width, color
%

% Find how many segments
xc = coast(:,1);
yc = coast(:,2);
L1 = isnan(xc);
nseg = sum(L1);

% Get start row index of each segment
i1 = find(L1)+1; % start

% Get end row of each segment
L2=L1;
L2(1)=0;
i2 = find(L2)-1;
i2 = [i2 ; length(xc)];

% Plot the figure
for n = 1:nseg;
   x = xc(i1(n):i2(n));
   y = yc(i1(n):i2(n));
   plot(x,y);
end





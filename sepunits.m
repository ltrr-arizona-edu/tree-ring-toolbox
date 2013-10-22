%  sepunits.m     Distance between 2 points on cartesian coord system.

%  Screen prompts for x, y coords of two points.
%  Distance is displayed.

x1=input('X-COORD OF FIRST POINT:  ');
y1=input('Y-COORD OF FIRST POINT:  ');
x2=input('X-COORD OF SECOND POINT:  ');
y2=input('Y-COORD OF SECOND POINT;  ');

delx=x2-x1;
dely=y2-y1;

dist=sqrt(delx^2  + dely^2)

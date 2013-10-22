function [d,c,ratio] =slp1(x1,y1,x2,y2)

% returns ratio of slopes of  y2 vs x2   to y1 vs x1

x1=[x1 ones(length(x1),1)];
c=x1\y1
pause



x2=[x2 ones(length(x2),1)];
d=x2\y2
pause

ratio=c(1)/d(1)
pause


pause





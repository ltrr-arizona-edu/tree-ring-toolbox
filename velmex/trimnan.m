function [x,yrx]=trimnan(x,yrx)
% trimnan:  Strip trailing and leading NaN
% [x,yrx]=trimnan(x,yrx);
x=trailnan(x);
mx=length(x);
yrx=yrx(1:mx);
x=flipud(x);
yrx=flipud(yrx);
x=trailnan(x);
mx=length(x);
yrx=yrx(1:mx);
x=flipud(x);
yrx=flipud(yrx);

function crnplot1(x,n,yrs,tit)
% crnplot1(x,n,yrs,tit)
% 
% Time series dual plot with (1) chronology and (2) sample size
% 
%
% D. Meko 2-29-96
%
%
%******************** IN ARGS ***********************************
%
% x (mx x 1)r time series of tree-ring index
% n (mn x 1)i time series of number of trees in sample
% yrs (1 x 2)i  first, last year of x and n
% tit (1 x 40)s   string info to place with gtext

close all

t = (yrs(1):yrs(2))';

a1=axes('Position',[0.13 0.45 0.775 0.415]);
h1=plot(t,x)
ylabel('Tree-Ring Index')
title(tit)
grid

v=axis;
a2=axes('Position',[0.13 0.11 0.775 0.25]);
stairs(t,n);

xlabel('Year')
ylabel('No. of Trees')
grid
set(gca,'YLim',[0 max(n)+2])

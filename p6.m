% p6.m   Potter's bivariate test on winter and summer ppt
%   Ft Bayard (suspect) vs Socorro (assumed good)

clg
subplot(111)

x1=pptw3(:,2);  % Bayard winter ppt, assumed to be good series.
y1= pptw2(:,2);  % Socorro winter ppt, assumed suspect
load tablepot;  % Table of critical points for T0

[T,I0,T0,dd,sig,T05]=potter(x1,y1,tablepot);


yrs=(pptw3(:,1));  % cv of years for analysis
n1=length(yrs);  % how many years passed to Potter.m
yrs(n1)=[];  % chop last year
yrkey=I0+yrs(1)-1;  % year of T0 = year before change in mean

plot(yrs,T,yrs,T05(ones(n1-1,1),:));
text(1960,T05,'95% sig point');
title('POTTERS BIVARIATE STATISTIC')
text( 1930,6,' WINTER PPT, SOCORRO VS BAYARD (GOOD)')
text(1940,5,['I0 = ',num2str(yrkey)]);
text(1940,4,['d = ',num2str(dd),' hundr in']);
xlabel('YEAR')
ylabel('T')

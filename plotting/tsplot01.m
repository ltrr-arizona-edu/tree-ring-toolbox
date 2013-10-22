function [h1,h2]=tsplot01(yr,x,n,tita,titb,ylaba,q,cshade)
% tsplot01: time series plot of an annual time series and sample size
% [h1,h2]=tsplot01(yr,x,n,tita,titb,ylaba,q,cshade);
% Last revised 5-9-99
%
% Time series is plotted in top axes, sample size below. User can shade
% y-axis range outside specified quantiles of x.  Horizontal line at median.
% Initial use: to display variation of annual precipitation averaged over all
% available stations in Santa Rita Experimental Range.  
%
%*** INPUT
%
% yr (? x 1)i  year vector for time series x
% x (? x 1)r time series
% n (? x 1)i sample size in each year (e.g., number of stations in reg-avg pcp
% tita{1 x 3} character strings of title info for top plot
%    {1} leftmost part of title (e.g., 'PPT')
%    {2} center part of title (e.g., 'Summer')
%    {3} rightmost part of title (e.g., 'July,August,September')
% titb{1 x 1} character string for title of sample size plot
%    note: this to be plotted in parens after "Sample Size" above lower axes
% ylaba{1 x 2} character strings of 2-part label for right axis of top plot
%      (e.g., {'PPT','(in)'};
% q (1 x 2)r lower and upper quantile defining shaded areas (e.g., [.1 .9])
% cshade(1 x 3)r color of shaded area (e.g., [.5 .5 .5])

close all;

%--- Set positions for plots
posb  = [.1 .1 .85 .22]; % for sample size
posa = [.1 .33 .85 .6]; % for time series

% Build axes
haxa = axes('position',posa);
haxb = axes('position',posb);

axes(haxa);
bar(yr,x);
ylabel([ylaba{1} '  ' ylaba{2}]);
set(gca,'XTickLabel',[],'Xgrid','on');
title([tita{1} '  ' tita{2} '  (' tita{3} ')']);

axes(haxb);
stairs(yr,n);
ylabel('Sample Size');
set(gca,'Xgrid','on',...
   'YLim',[max([0,min(n)-range(n)/10])   max(n)+range(n)/10]);
legend(titb{1});



disp('this')






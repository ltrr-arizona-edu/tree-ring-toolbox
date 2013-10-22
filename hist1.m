% hist1.m  make a data file that grapher can use to plot a histogram.
%     Also gives matlab screen plot of the histogram.
%
%  dedicated version 5-8-91, where X is ARFIT.DAT and j=3, and
%   s=0:10:70
%
%***************   PRELOADS  ********************************

% X .... the data array, one column of which is the variable to be treated

%*******************  SCREEN PROMPTS  *******************

% s .... row vector of class intervals (eg: 0:10:70)

% j ...  column of D holding the data for histogram

%*****************************************************************

clg


s=input('VECTOR OF CLASS INTERVALS: EG: 0:10:70  ? ')
j=input('COLUMN OF INPUT DATA ARRAY TO BE HISTOGRAMMED? ')

s=0:10:70;

%********** Make matlab screen plot of histogram

hist(X(:,j),s);
xlabel('PERCENT VARIANCE EXPLAINED')
ylabel('NUMBER OF SITES')
title('MATLAB HISTOGRAM')
pause

[N,x]=hist(X(:,j),s);   % Store histogram columns
[XX,YY]=bar(x,N);    %  Convert to bar-plot variables

plot(XX,YY);
ylabel('FREQUENCY');
xlabel('PERCENT VARIANCE EXPLAINED');
title('HISTOGRAM USING PLOT(XX,YY)');

%***************  Store the two variables for later grapher histogram plot

data=[XX YY];

save hist.dat data /ascii

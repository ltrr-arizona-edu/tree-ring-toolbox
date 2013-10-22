% p2.m  project question  2:  scatterplot and time series

% Assumes have all tree-ring series in a matrix x
% the predictand in a col vector y

years=[1929 1979];  % cv of years for plot
yrtrees=(1700:1979)';  % cv of years for full tree-ring array
yrpred = (1929:1989)';  % cv of years for predictand array

L1 = yrtrees >= years(1) & yrtrees <= years(2);
L2 = yrpred >= years(1) & yrpred <= years(2);

x1=x(L1,1);  % desired tree-ring series for plot.  
y1=y(L2,1);  % desired subyears of streamflow
yr=(years(1):years(2))';  % cv of years for plot

[b,r,yhat]=lintrnd(y1,x1);  % get linear regression line and r

plot(y1,x1,'*',y1,yhat);
title('MIMBRES TRI VS GILA R FLOW NEAR GILA, NM')
xlabel('WY FLOW (ACRE-FT))')
ylabel('TREE-RING INDEX');
text(200000,.6,['r = ',num2str(r)]);
pause

y1=(y1-mean(y1))/std(y1);
x1=(x1-mean(x1))/ std(x1);

plot (yr,x1,'-',yr,y1,'--');
title ('MIMBRES TRI (SOLID) AND WY FLOW (DASHED)')
xlabel('YEAR')
ylabel('Z-SCORES');


% p3prep.m   prepares arrays of seasonalized ppt and temp

% Assumes lunp,socp,bayp, lunt,soct,bayt are the monthly
% ppt and temp series for 3 stations, all covering same period,
% and in workspace

% Make arrays with col 1 as year, and cols 2-4 as data for
% Luna, Socorro, and Bayard, in that order.  Call the 
% arrays PPTW, PPTS,  TMPW,  TMPS -- ppt and temperature
% for winter and summer.  Let Nov-Apr be winter, June-Sept be
% summer.


pptw1 = seaspt(lunp,11,4,[1905 1988],1);
pptw2 = seaspt(socp,11,4,[1905 1988],1);
pptw3 = seaspt(bayp,11,4,[1905 1988],1);
tmpw1 = seaspt(lunt,11,4,[1905 1988],2);
tmpw2 = seaspt(soct,11,4,[1905 1988],2);
tmpw3 = seaspt(bayt,11,4,[1905 1988],2);


ppts1 = seaspt(lunp,6,9,[1905 1988],1);
ppts2 = seaspt(socp,6,9,[1905 1988],1);
ppts3 = seaspt(bayp,6,9,[1905 1988],1);
tmps1 = seaspt(lunt,6,9,[1905 1988],2);
tmps2 = seaspt(soct,6,9,[1905 1988],2);
tmps3 = seaspt(bayt,6,9,[1905 1988],2);



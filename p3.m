% p3a.m  timeseries plots of ppt for stations lun, soc, bay on 1 plot

% assumes that have used p3prep.m to make the seasonalized climate series

%***** Plots for winter precipiation

plot(pptw1(:,1),pptw1(:,2),'-',  pptw2(:,1),pptw2(:,2),'--',...
     pptw3(:,1),pptw3(:,2),':');
title('WINTER (NOV-APR) PRECIPITATION AT THREE STATIONS')
text(1930,1600,'solid = Luna RS')
text(1930,1400,'dashed = Socorro')
text(1930,1200,'dotted = Bayard')
xlabel('YEAR')
ylabel('PPT (HUNDREDTHS OF INCHES)')
pause

rwp=corrcoef([pptw1(:,2)  pptw2(:,2)   pptw3(:,2)]) ; % correlation matrix
%   for winter precipitation at the three stations.


%****** Similar plots for winter temperature


plot(tmpw1(:,1),tmpw1(:,2),'-',  tmpw2(:,1),tmpw2(:,2),'--',...
     tmpw3(:,1),tmpw3(:,2),':');
title('WINTER (NOV-APR) AVERAGE TEMPERATURE AT THREE STATIONS')
text(1960,490,'solid = Luna RS')
text(1960,480,'dashed = Socorro')
text(1960,470,'dotted = Bayard')
xlabel('YEAR')
ylabel('TEMPERATURE (TENTHS OF DEGREES F)')
pause

rwt=corrcoef([tmpw1(:,2)  tmpw2(:,2)   tmpw3(:,2)]) ; % correlation matrix
%   for winter temperature at the three stations.



%***** Similar plots for summer precipiation

plot(ppts1(:,1),ppts1(:,2),'-',  ppts2(:,1),ppts2(:,2),'--',...
     ppts3(:,1),ppts3(:,2),':');
title('SUMMER (JUN-SEPT) PRECIPITATION AT THREE STATIONS')
text(1930,1600,'solid = Luna RS')
text(1930,1400,'dashed = Socorro')
text(1930,1200,'dotted = Bayard')
xlabel('YEAR')
ylabel('PPT (HUNDREDTHS OF INCHES)')
pause

rsp=corrcoef([ppts1(:,2)  ppts2(:,2)   ppts3(:,2)]) ; % correlation matrix
%   for summer precipitation at the three stations.



%****** Similar plots for summer temperature


plot(tmps1(:,1),tmps1(:,2),'-',  tmps2(:,1),tmps2(:,2),'--',...
     tmps3(:,1),tmps3(:,2),':');
title('SUMMER (JUNE-SEPT) AVERAGE TEMPERATURE AT THREE STATIONS')
text(1960,780,'solid = Luna RS')
text(1960,770,'dashed = Socorro')
text(1960,760,'dotted = Bayard')
xlabel('YEAR')
ylabel('TEMPERATURE (TENTHS OF DEGREES F)')
pause

rst=corrcoef([tmps1(:,2)  tmps2(:,2)   tmps3(:,2)]) ; % correlation matrix
%   for summer temperature at the three stations.


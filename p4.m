% p4.m  scatterplot of  seasonalized p vs t at Ft Bayard

% winter first

%  Call for regression line

[b,r,yhat]=lintrnd(pptw3(:,2),tmpw3(:,2));

plot(pptw3(:,2),tmpw3(:,2),'*',pptw3(:,2),yhat);
title('FT BAYARD: NOV-APR TEMPERATURE VS PPT')
xlabel('TEMPERATURE (TENTHS OF DEGREE F)')
ylabel('PPT (HUNDREDTHS OF INCHES)')
text(800,480,['r = ',num2str(r)]);

pause


% then summer

[b,r,yhat]=lintrnd(ppts3(:,2),tmps3(:,2));

plot(ppts3(:,2),tmps3(:,2),'*',ppts3(:,2),yhat)
title('FT BAYARD: JUNE-SEPT TEMPERATURE VS PPT')
xlabel('TEMPERATURE (TENTHS OF DEGREE F)')
ylabel('PPT (HUNDREDTHS OF INCHES)')
text(1200,760,['r = ',num2str(r)]);
pause

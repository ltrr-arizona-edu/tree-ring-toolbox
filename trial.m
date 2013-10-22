clg
x=A(:,2);
y=A(:,3);
yr=A(:,1);
subplot(211);
plot(yr,x,'-',yr,y,'--');
xlabel('year')
ylabel('ppt (in)');
pause

[b,r,yhat]=lintrnd(x,y);

plot(x,y,'*',x,yhat,'-')
title('SCATTERPLOT');
text(.5,.5,['r = ',num2str(r)]);
pause
subplot(111);

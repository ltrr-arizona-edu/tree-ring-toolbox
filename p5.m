% p5.m   double mass plots for 1905-1989 winter and summer ppt,
%  at stations 1, 2 3:  Luna, Socorro, Bayard

clg

% Begin with winter analysis, then summer

subplot(111)

% Luna vs Soccoro

x1=pptw2(:,2);
y1=pptw1(:,2);
x2=cumsum(x1);
y2=cumsum(y1);
[tau,crit]=mannkend(x1,y1);
subplot(221)
plot(x2,y2,'*')
ylabel('LUNA WINTER PPT')
xlabel('SOCORRO WINTER PPT')

xp1=min(x2)+(max(x2)-min(x2))/4; % plotting x point for tau
yp1=min(y2)+3*(max(y2)-min(y2))/4; % plotting x point for tau
text(xp1,yp1,['tau=',num2str(tau)]);
yp2=min(y2)+3*(max(y2)-min(y2))/5; % plotting x point for p05
text(xp1,yp2,['p05=',num2str(crit)]);


% Luna vs Bayard


x1=pptw3(:,2);
y1=pptw1(:,2);
x2=cumsum(x1);
y2=cumsum(y1);
[tau,crit]=mannkend(x1,y1);
subplot(222)
plot(x2,y2,'*')
ylabel('LUNA WINTER PPT')
xlabel('FT BAYARD WINTER PPT')

xp1=min(x2)+(max(x2)-min(x2))/4; % plotting x point for tau
yp1=min(y2)+3*(max(y2)-min(y2))/4; % plotting x point for tau
text(xp1,yp1,['tau=',num2str(tau)]);
yp2=min(y2)+3*(max(y2)-min(y2))/5; % plotting x point for p05
text(xp1,yp2,['p05=',num2str(crit)]);



% Socorro vs Bayard


x1=pptw3(:,2);
y1=pptw2(:,2);
x2=cumsum(x1);
y2=cumsum(y1);
[tau,crit]=mannkend(x1,y1);
subplot(223)
plot(x2,y2,'*')
ylabel('SOCORRO WINTER PPT')
xlabel('FT BAYARD WINTER PPT')

xp1=min(x2)+(max(x2)-min(x2))/4; % plotting x point for tau
yp1=min(y2)+3*(max(y2)-min(y2))/4; % plotting x point for tau
text(xp1,yp1,['tau=',num2str(tau)]);
yp2=min(y2)+3*(max(y2)-min(y2))/5; % plotting x point for p05
text(xp1,yp2,['p05=',num2str(crit)]);
pause

% Repeat the above, but for summer ppt


subplot(111)

% Luna vs Soccoro

x1=ppts2(:,2);
y1=ppts1(:,2);
x2=cumsum(x1);
y2=cumsum(y1);
[tau,crit]=mannkend(x1,y1);
subplot(221)
plot(x2,y2,'*')
ylabel('LUNA SUMMER PPT')
xlabel('SOCORRO SUMMER PPT')

xp1=min(x2)+(max(x2)-min(x2))/4; % plotting x point for tau
yp1=min(y2)+3*(max(y2)-min(y2))/4; % plotting x point for tau
text(xp1,yp1,['tau=',num2str(tau)]);
yp2=min(y2)+3*(max(y2)-min(y2))/5; % plotting x point for p05
text(xp1,yp2,['p05=',num2str(crit)]);


% Luna vs Bayard


x1=ppts3(:,2);
y1=ppts1(:,2);
x2=cumsum(x1);
y2=cumsum(y1);
[tau,crit]=mannkend(x1,y1);
subplot(222)
plot(x2,y2,'*')
ylabel('LUNA SUMMER PPT')
xlabel('FT BAYARD SUMMER PPT')

xp1=min(x2)+(max(x2)-min(x2))/4; % plotting x point for tau
yp1=min(y2)+3*(max(y2)-min(y2))/4; % plotting x point for tau
text(xp1,yp1,['tau=',num2str(tau)]);
yp2=min(y2)+3*(max(y2)-min(y2))/5; % plotting x point for p05
text(xp1,yp2,['p05=',num2str(crit)]);



% Socorro vs Bayard


x1=ppts3(:,2);
y1=ppts2(:,2);
x2=cumsum(x1);
y2=cumsum(y1);
[tau,crit]=mannkend(x1,y1);
subplot(223)
plot(x2,y2,'*')
ylabel('SOCORRO SUMMER PPT')
xlabel('FT BAYARD SUMMER PPT')

xp1=min(x2)+(max(x2)-min(x2))/4; % plotting x point for tau
yp1=min(y2)+3*(max(y2)-min(y2))/4; % plotting x point for tau
text(xp1,yp1,['tau=',num2str(tau)]);
yp2=min(y2)+3*(max(y2)-min(y2))/5; % plotting x point for p05
text(xp1,yp2,['p05=',num2str(crit)]);


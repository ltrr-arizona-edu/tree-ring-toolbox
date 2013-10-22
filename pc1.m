clg;
x1=[6 10 9 5 10]';
x2=[82 78 78 82 80]';
x3=[.2 1.8 1 0.5 1.5]';
year=[1 2 3 4 5]';
subplot(221);
axis([0 6 0 12]);
plot(year,x1,year,x2-80,year,x3);
axis([0 12 0 3]);
title('X1, X2, X3 vs Year');
xlabel('Year'), ylabel('Value');
subplot (222);
plot(x1,x3,'x');
title('X3 vs X1');
xlabel('PPT'), ylabel('Ring');
subplot (223);
axis([0 12 75 85]);
plot(x1,x2,'x')
title('X2 vs X1');
xlabel('PPT'), ylabel('Temp');
subplot (224);
axis([75 85 0 3]);
plot(x2,x3,'x');
title('X3 vs X2');
xlabel('Temp'), ylabel('Ring');

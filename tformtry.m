function tformtry(x1,y1);

m1=length(x1);
m2=length(y1);
%yrsx=input('START, END YEAR OF TREE-RING ARRAY: ');
%yrsy=input('START, END YEAR OF PREDICTAND ARRAY: ');

yrsx=[1700 1979];
yrsy=[1929 1989];


% Find overlap period

yrsover = [max([yrsx(1) yrsy(1)])  min([yrsx(2)  yrsy(2)])];
yrx = (yrsx(1):yrsx(2))';
yry = (yrsy(1):yrsy(2))';

Lx= yrx >= yrsover(1) & yrx <= yrsover(2);
Ly= yry >= yrsover(1) & yry <= yrsover(2);

x2=x1(Lx);  y2=y1(Ly);

plot(x2,y2,'*')
title('GILA ANNUAL FLOW VS MIMBRES TRI')
xlabel('INDEX');
ylabel('ACRE-FT');
pause

x3=x2;
y3=log(y2);
plot(x3,y3,'*')
title('LOG GILA ANNUAL FLOW VS MIMBRES TRI')
xlabel('INDEX');
ylabel('LOG ACRE-FT');
pause


x3=x2;
y3=1 ./ y2;
plot(x3,y3,'*')
title('RECIPROCAL GILA ANNUAL FLOW VS MIMBRES TRI')
xlabel('INDEX');
ylabel('1 / ACRE-FT');
pause


x3= 4 .^x2;
y3=y2;
plot(x3,y3,'*',exp(x2),y3,'+')
title('GILA ANNUAL FLOW VS SQUARED MIMBRES TRI')
xlabel('SQUARED INDEX');
ylabel('ACRE-FT');
pause'


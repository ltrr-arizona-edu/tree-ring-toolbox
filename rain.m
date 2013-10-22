% ppt1.m   makes a 4-season bar chart of a monthly ppt series

%****************   PRELOADS   *******************************************
%
% P		Monthly ppt.  First col a year, then 13 data values
%
% A      Seasons mask.  ?  x  24. Tells which months go into each 
%		season.  Element 1 is jan of year t-1, element 24 is dec
%		of year t.  For example, here is a mask for an Sept-Aug
%		total:
%
%			0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0
%
%**************  OTHER VARIABLES   ******************************************
% 
% m1,n1	total # rows and columns in input P	
% m2,n2	total # of seasonal series to make; and a moot 24
%
% P1.....duped and lagged 24-col version of P, corrected for integer
%		input format
% P2.....seasonal and annual total series
% P3.....cumulative total seasonal series (for stacked bars)
% yr.....cv, "tree-ring" year
%
%
%*****************  COMPUTE SEASONAL TOTAL SERIES  *************************

[m1,n1]=size(P);
[m2,n2]=size(A);

yr = P(2:m1,1);

P1=[P(1:m1-1,2:13) P(2:m1,2:13)];  % make 24-col monthly ppt array

P1= 0.01 * P1;   % correct for integer format of trl input ppt

P2 = P1 * A';   % compute the seasonal-total series


%****************  COMPUTE STACKED BAR SERIES ************************

P3 = zeros(m1-1,m2-1);  %  pre-allocate
P3(:,1)=P2(:,1);  % First season stands alone

for i = 2:m2-1
	P3(:,i)=P3(:,i-1) + P2(:,i);  % form fall+winter, fall+win+spr, etc
end


if any ((P2(:,m2) - P3(:,m2-1)) > .00001);  % Flag if annual total wrong
	disp('Annual Totals Dont Jive');
	pause
end

%**************** MAKE A STACKED BAR ********************************
%

[x1,y1]=stairs(yr,P2(:,1));
[x2,y2]=stairs(yr,P2(:,2));
[x3,y3]=stairs(yr,P2(:,3));
[x4,y4]=stairs(yr,P2(:,4));
[x5,y5]=stairs(yr,P2(:,5));
m3=length(x1);

y1(1)=y1(2);
y1(m3)=y1(m3-1);
y2(1)=y2(2);
y2(m3)=y2(m3-1);
y3(1)=y3(2);
y3(m3)=y3(m3-1);
y4(1)=y4(2);
y4(m3)=y4(m3-1);
y5(1)=y5(2);
y5(m3)=y5(m3-1);

clg

axis([1850 2000 0 15]);

%plot(x1,y1,x2,y2,x3,y3,x4,y4);
%title('TUCSON U OF A');
%xlabel('YEAR');
%ylabel('INCHES');
%pause

subplot(221);
plot(x1,y1);
title('SEP-NOV');
grid

subplot(222);
plot(x2,y2);
title('DEC-FEB');
grid

subplot(223);
plot(x3,y3);
title('MAR-MAY');
grid

subplot(224);
plot(x4,y4);
title('JUN-AUG');
grid

pause;

axis([1900 2000 0 25]);

subplot
plot(x5,y5);
grid
xlabel ('YEAR')
ylabel ('PRECIPITATION (INCHES)')

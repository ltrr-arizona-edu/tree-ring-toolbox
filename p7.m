%  p7.m  tables of basic statistics for selected variables

%  Put vbl number in col 1, means in col2, stdev in col 3, 
%  skew in col 4,  first-order r in col 5
% Use common period 1929-79 covered by all data

% variable definitions

% 1 winter ppt station 1
% 2 winter ppt station 2
% 3 winter ppt station 3

% 4 winter temp station 1
% 5 winter temp station 2
% 6 winter temp station 3

% 7 summer ppt station 1
% 8 summer ppt station 2
% 9 summer ppt station 3

% 10 summer temp station 1
% 11 summer temp station 2
% 12 summer temp station 3

% 13 predictand 
% 14-18  tree-ring index 1929-79



%****** Put all series in one array
%****** Start with the winter ppt and temp series
%     A scaling factor is included because ppt and temp unit are
%     stored in hundredths of inches and tenths of degrees.

yr1=1929;
yr2=1979;


L1=pptw1(:,1) >= yr1 & pptw1(:,1) <= yr2;

X1(:,1:6) = [pptw1(L1,2)*.01  pptw2(L1,2)*.01  pptw3(L1,2)*.01  ...
             tmpw1(L1,2)* 0.1  tmpw2(L1,2)*0.1  tmpw3(L1,2)*0.1];

% Next build columns of summer ppt and temp series

L1=ppts1(:,1) >= yr1 & ppts1(:,1) <= yr2;

X1(:,7:12) = [ppts1(L1,2)*.01  ppts2(L1,2)*.01  ppts3(L1,2)*.01  ...
             tmps1(L1,2)* 0.1  tmps2(L1,2)*0.1  tmps3(L1,2)*0.1];

%*****  Get appropriate years of streamflow series
%  This is stored as variable y, a cv covering year 1929-1989

yy=(1929:1989)';  % make cv of years of y
L1=yy >= yr1 & yy <= yr2;  %  pointer to desired years of y
X1(:,13)=y(L1);    % Subset of y, covering years yr1 to yr2

%***** Get appropriate years of the five tree-ring series
%*****  These series are stored in set2.mat as x

yy=(1700:1979)';   % cv of years in x
L1 = yy >= yr1 & yy <= yr2;   % pointer to desired years of x
X1(:,14:18) = x(L1,:);


%********      Compute statistics

% functions mean.m, std.m and skew.m all operate on matrices, while
% acf.m will require a loop

[m1,n1]=size(X1);

for i=1:n1
	[r pr]=acf(X1(:,i), 10);
	r1(i)=r(2);
end

[g,pr]=skew(X1);

SX = [(mean(X1))'  (std(X1))'  g'  r1']

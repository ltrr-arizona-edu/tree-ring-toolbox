% run2.m     compute and plot # of cells drier than a threshold

clear d yr m n m1 m2 n1 n2 X1 X2 x3 c xx yy
clear i i1 i2 f


d = -0.6;   % threshold number of standard deviations
yr = X(:,1);

[m,n]=size(X);
n1=n-1;    % number of cells
m1=m;      % number of years

[m3,n3]=size(V);

X1=X(:,2:n);  	% year col from X not included
%a = mean(X1); % rv of means of cell-ave indices
%b = std(X1);  % rv of stand devs

[Y,I]=sort(X1);  % sort each cell-ave tri from driest to wettest year

XS=[Y(fix(m/3),:);  Y(ceil(m/3),:)];
a=mean(XS);   % tercile cutoff point for each index series

%X1=X1 - a(ones(m1,1),:);  % subtract means from columns
%X1=X1 ./  b(ones(m1,1),:);  % complete conversion to z scores

% Compute array with entries 1 if cell drier than trheshld, 0 otherwise

X2=X1 < a(ones(m1,1),:);

%X2 = X1 - d(ones(m1,1),ones(n1,1)) <= 0;
x3=sum(X2');  %  rv of # of cells drier than threshld each year
x3=x3';  % convert to cv

% Compute median # of cells below thhold.  Then make stair plot of
% # of regions below threshold, with median # as base line.

c=-9;
%c = -1.0 * median(x3);   % negate to make dry side on bottom of line
clg,  subplot;
subplot(211);
for i = 1:4;  %m3
	clear xx yy	
	i1 = IY(i,1) - IYRS(1) + 1;
	i2 = IY(i,2) - IYRS(1) + 1;
	[xx,yy]= stairs(yr(i1:i2),x3(i1:i2));
	m2 = length(xx)+1;
	yy(1) = yy(2);
	yy(m2) = yy(m2-1);
   xx=xx-0.5;
	xx(m2) = xx(m2-1) + 0.5;

	f = -1.0;
	yy = yy .* f(ones(m2,1),:);
	axis(V(i,:));
	
	plot(xx,yy,'-r',xx,c(ones(m2,1),:),'-g');
	xlabel('YEAR');
	ylabel('- NUMBER OF CELLS');

%  Save files for plotting with grapher

%	savedat=[xx yy c(ones(m2,1),:)];
%	eval(['save C:\projs\ac1\file',int2str(i),'.dat savedat -ascii']);


	%meta c:\work\4f

	pause;
end

%  save cell counts in arrays that can be imported to spreadsheet

x4=x3-abs(c);
ll=x4<=0;
x4(ll)= zeros(sum(ll),1);


%cnt=[yr(1:100) x4(1:100) yr(101:200) x4(101:200) yr(201:300) x4(201:300)];
%cnt1=[yr(301:363)  x4(301:363)];

%save c:\projs\ac1\cnt.dat cnt -ascii
%save c:\projs\ac1\cnt1.dat cnt1 -ascii




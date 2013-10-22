function M=findmiss(X);

%  Finds years and months with years of missing data, as
%  marked by values exceeding hard-coded  "misscode".

%  Operates on monthly precipitation or temperature arrays.
%  Requires that no rows (complete years) be missing from the arrays.

% Makes a 2-col matrix:

% row 1, col 1: numeric station id
% row 1, col 2: data type (1=ppt, 2=temp)
% row 2:  first and last year of data matrix

% Remaining rows:

% col 1:  month with data missing
% col 2:  year with data missing


id=input('8-DIGIT STATION ID:  ');
itype=input('1 for PPT, or 2 for TEMP:   ');

[m1,n1]=size(X);
I=  X>=9998;
yrgo=X(1,1);
yrstop=X(m1,1);
M=[id itype; yrgo yrstop];

for j=1:12
	I1=I(:,j+1);
	irow=find(I1==1);
	yr=irow+yrgo-1;
	M=[M; j(ones(length(yr),1),:) yr];
end


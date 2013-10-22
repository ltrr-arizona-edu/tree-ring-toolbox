% lyrs.m  write to a file  specified years of data from two arrays
%  limited to one variable from each array

% May want to edit write filename, header, and col headings
% Writes: header line, blank line, col headings, blank line, and 
%     table of data


%******************  preloads *********************
% k ... cv of years to be listed
% w ... col from X to be used
% v ... col from Y to be used
% b ... beginning year of array X
% c ... beg year of Y
% X ... data array, possibly with year as first col
% Y ... data array
%
%*******************************************************

r=0;
b1=b-1;
c1=c-1;
k1=k-b1(ones(length(k),1),:);
k2=k-c1(ones(length(k),1),:);

clc

x1=X(k1,w);
y1=Y(k2,v);

clc
home
fprintf('c:\scratch\hey.dat','20-yr Values Keyed on Colo R. Dry Years\n\n')
fprintf('c:\scratch\hey.dat','year    Colo   Four Rivs\n\n');

for i=1:length(k)
	fprintf('c:\scratch\hey.dat','%5.0f %8.3f % 8.3f\n',k(i),x1(i),y1(i))
end


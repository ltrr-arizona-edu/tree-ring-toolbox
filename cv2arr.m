% cv2arr.m     column vector to an array

% D. Meko,    9-26-91

% Converts several time series strung out in a single column vector
% to an array in which the time series are separate columns.

% Assumes input column vector contains two or more time series, each
% covering exactly the same  years.


yr1 = input('FIRST YEAR OF DATA:  ');
yr2 = input('LAST YEAR OF DATA:  ');

n2 = input('NUMBER OF TIME SERIES:  ');

m2=yr2-yr1+1;  % number of years


Y=zeros(m2,n2);  % preallocate target array.

%******  LOAD THE STRUNG-OUT VECTOR, IF NOT ALREADY LOADED  *******

xin=input('HAS THE SOURCE VECTOR BEEN LOADED? Y/N [N] ','s');
if isempty(xin), xin='N'; end;
if xin=='N'
	disp('LOAD THE VECTOR winter.dat')
	keyboard
else
end

Y(:)=summer;

function Yhat = replacem(Y,Z,U)

% Replace monthly ppt or temp data.

% Given a ppt series or temp series, a same-sized array of estimated data,
% and a pointer array to years and months, replaces values of the ppt or
% temp array with values from the estimated data.

% Assumes that you have previously run estppt.m or esttmp.m to get
% the arrays Z and U.

%***********  INPUT ARGS

% Y (? x 13) the ppt or temperature array
% Z (? x 13) the corresponding array containing values estimated from 
%		estppt.m or esttmp.m
% U (? x 5)  a pointer array telling which years and months of Y are
%	to be replaced with values from Z

%********   OUTPUT ARGS

% Yhat (? x 13)  same as Y, except with selected values replaced from Z


[m1,n1] = size(U);  
num=m1-2;  % Number of monthly values to replace.  Note that
%	first two rows of U contain bookkeeping info.

T = U (3:m1,:);  % U, not including the two header rows
Yhat=Y;  % Initialize 

for i = 1:num;   %  Loop for each value to be replaced
	row = T(i,2) - Y(1,1) + 1;  % convert year to row pointer
	col = T (i,1) +1;  % Convert month to col pointer
	Yhat (row,col) = Z(row,col);
end

% shap2.m  Unshapes a matlab matrix into standard time series form.
% 
% The new matrix may have multiple rows per year, and will have
% dummy variables to make the number of variables evenly divisible
% by the number of variables in a record.

% X2 is the original matlab matrix. If too many cols to fit in a standard
%	 255 char record, will not be able to simply save as ascii for later 
% 	 use in Fortran.

% X1 is the new, reshaped matlab matrix.  Will not include the year.


l = input('Number of desired rows per year in target array: ')
I1 = input('First year: ')
I2 = input('Last year: ')

X2(:,1)=[];  % strip the year column off X2

m=I2-I1+1  % number of years, which ==  number of rows in X2

[m,n] = size(X2);

A1= zeros(n/l,m*l);

B1=X2';
A1(:)=B1;
X1=A1';

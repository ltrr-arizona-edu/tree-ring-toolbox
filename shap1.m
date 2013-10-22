% shap1.m  Shapes a matrix with approp rows as years and cols as vars
% Assumes you have run SHAP1.F to pre-treat array, then
%  loaded array into matlab and renamed as X1.  
% Gives X2.DAT as output array.
% Gives IYEAR as year col vector

n = input('Number of desired columns, or variables: ')
I1 = input('First year: ')
I2 = input('Last year: ')



m=I2-I1+1  % number of years, which == desired number of rows

X2=zeros(m,n);

IYEAR = [I1:I2]';

% assume the input array that needs to be rearranged is X1
% assume you have pre-loaded X1
% assume the created output array in X2


A1=X1';

B1 = ones(n,m);

B1(:)=A1;

X2=B1';

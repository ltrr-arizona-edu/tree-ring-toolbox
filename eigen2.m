function [chi,df,T] = eigen2(x,L)

% Chi squared test for zero intercorrelation;  Table of eigenvector
% pct variance explained.

% Feb 29, 1992, by D. Meko

%*********   INPUT ARGS
%
% x  (m1 x n1) data array,  m1 rows, n1 cols
% L  (diagonal array of ordered eigenvalues of x
%     Obtained previously using eigen1.m

%*******   OUTPUT ARGS
%
% chi (1 x 2)   test statistic and deg freedom  for 
% 	 Bartlett's sphericity test (Cooley and Lohnes 1971, p. 103)
%   H0: variables are already uncorrelated before pca
%   Evaluate from chi squared table.
% df (1 x 1)  degrees of freedom for chi squared test on chi
% T (n1 x 4)  table of following quantities for each chronology (row)
%
%	col 1 --  component number
%  col 2 --  eigenvalue (from diagonal of L)
%  col 3 --  pct variance accounted for by this component
%  col 4 --  cumulative pct variance accounted for by components
%            1 through this component. 
%            

%**************

t = diag(L);     % eigenvalues
ts=sum(t) ;  % pct variance expld by eigenvectors
[N,p]=size(x);  

T=zeros(p,4);  % Initialize Table

DR = det(corrcoef(x));

T(:,1:4) = [(1:p)'  t  100*t/ts  100*cumsum(t/ts)];


%**********  Begin sphericity test calcs
%
% First, for H0 that det of population correlation matrix is
% already 1 (i.e., variables are already uncorrelated without
% converting to components

chi = - ((N-1) - (1/6)* (2*p+5)) *  log(DR); % test criterion
%		df = 0.5 (p-squared -p)
df = 0.5 * (p^2-p);


function [x,p]=chi1(C)
% chi-squared and corresp limiting values at 5%, 1%  levels

%**************************  INPUT ARGS
%
% C (mC x mC) -  contingency table


%********  OUTPUT ARGS
%
% x (1 x 1)  -  chi-squared statistic
% p (1 x 2) - corresp table chi-sq at 5%, 1% levels


%*******************   OTHER USER-WRITTEN FUNCTIONS NEEDED
%
% chi_5_1.mat  - chisquared table, 5% and 1% levels; must be global
%		in workspace


%************  USE EXAMPLE
%
% Null hypoth: rw variations between core a and core b are 
% no different than expected with no relationship (no crossdating).
% If x> p(1), reject at 5%,  if x>p(2), reject at 1%


%********************** GO

global chi_5_1  % function must access table of chi-sq values


[m,n]=size(C);
if (m ~= n), error('C not square'), end;

% Say sample groups are cols and master groups are rows
% Compute marginal totals
S1=sum(C);  % marginal totals for sample groups
M1=(sum(C'))'; % marginal totals for master groups

% Compute expected number in each cell
E = (M1 * S1) / (sum(M1));

% Compute statistic
xcell=((abs(C-E)) .^2) ./ E;
x=sum(sum(xcell));

% Compute degrees of freedom
df=(m-1)* (n-1);


% Get 5%, 1% levels of chi-sq for this no. df
p=table1(chi_5_1,df)

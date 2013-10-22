function [A,k1,vrat,C,S2]= arwhite2(X,nhi,k2);
% arwhite2: fits cols of a matrix to ar models, culling non-NaN segment
% CALL: [A,k1,vrat,c,s]= arwhite2(X,nhi,k2);
%
% Meko 10-2-97
%
%********************** IN *********************
%
% X (mX x nX)r  time serie matrix.  Each series may cover diff periods
% nhi (1 x 1)i  order of ar model. Either the highest 
%   order to consider (if k2==1) or the order to fit (if k2==2)
% k2 (1 x 1)i option for order to fit (see nhi)
%
%********************  OUT ***********************
%
% A (mA x nA)r  ar residuals matrix, same size of X
% k1 (nA x 1)i  order of ar model fit (might be 0)
% vrat (nA x 1)r ratio of variance of resids to variance of orig
%
%****************** NOTES ************************
%
% Modeled on whit1.m.  arshite2.m calls whit1.m for each time series
% in matrix.  
%
% Handles X matrices that has NaNs for various elements of cols of X
%
% No year col in X or A.  I assume the caller has a year vector
% in his pocket on the side
%
% User may want to truncate beginning years of chrons before modeling
% if accectpable sample size (ass) too small.  see truncass.m


% Size X and A
[mX,nX]=size(X);
t = (1:mX)';  % dummy 'year' column

anan = NaN;
A = anan(ones(mX,1),ones(nX,1)); % to hold residuals
k1 = anan(ones(nX,1),:); % to hold model order
vrat = anan(ones(nX,1),:);
C = anan(ones(nX,1),ones(6,1)); % to hold ar coefs
S2 = C; % to hold 2-standard errors of ar coefs

% Check k2
if k2 ~=1  & k2~2;
   error('k2 must be 1 or 2');
end


L = ~isnan(X); % marks valid data

clc
disp('MODELING SERIES');
for n = 1:nX;
   disp(['Series # ' int2str(n)]);
   L1 = L(:,n);
   x = X(L1,n);
   [e,kk1,vvrat1,arcs] = whit1(x,nhi,k2);
   A(L1,n) = e;
   vrat(n)=vvrat1;
   k1(n)=kk1;
   if kk1==0;
      % no action , leave row for C and S2 as NaN
   else
      C(n,1:kk1)=arcs(1,:);
      S2(n,1:kk1) = arcs(2,:);
   end
end

           
   

function [BC,D]= pdsicnt3(B,b,mask)
% Tabulates number of gridpoints out of 153 with reconstructed pdsi
% below a specified threshold.  For Ed Cooks gridpoint reconstruction.
% Written with pdsicnt2.m as template.
%
% D Meko 11-2-94
%
%*******   INPUT ARGS
%
% B (mB x nB)  reconstructed PDSI, mB years and nB-1 gridpoints; col 1=yr
% b(1 x nB-1) threshold drought value for each of nB-1 gridpoints.  
% 	 This would typically be a constant vector of some value.  Say 
% 	 b(1:nB-1)=-3.0.  Then only years with PDSI < -3.0 will classify as
%   drought years.
% mask (1 x nB-1):  1 if include this site, 0 if skip

%**********  OUTPUT ARGS   **************

% BC (mB x 1)   recon PDSI count of number of points in drought each year
% D (1 x nB-1)   summary information on drought for each gridpoint
%   row 1:  number of years classified in drought- recon PDSI


%************  SIZE SOME ARRAYS
[mB,nB]=size(B);


%************  PREALLOCATE
BC=zeros(mB,1);
yrB = B(:,1);
B(:,1)=[];  % cut year col off B
[mB,nB]=size(B);
D=zeros(1,nB);

%*****************  EXPAND THRESHOLD DROUGHT VECTOR INTO A MATRIX
BT = b(ones(mB,1),:);



%**************** Count number of years in drought
LB1=B<=BT;  
D = sum(LB1);
BC =  (sum(LB1'))';



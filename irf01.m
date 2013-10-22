function [c,lag,s]=irf01(Z,m,na)
% irf01: sign of largest of first m irf weights, and whether sig at 99%
% CALL: irf01
%
% Sign of largest of first m impulse response weights , and 
% % Lags 0-4
%
% Meko 10-30-97
%
%******************  IN ************
%
% Z (mZ x 2)r  output (col 1)  and input (col 2) time series
% m (1 x 1)i  compute to this many lags
% na (1 x 1)i  ar order to whiten input and output before computing irf
%
%**************  OUT **************
%
% c (1 x 1)r  largest irf weight
% lag (1 x 1)i  lag of that weight
% s (1 x 1)L   significant (1) or not (0) at 99%


% Compute impulse response function
plt=0;
[IR,R,CL] = cra(Z,m,na,plt) ;

% Find largest absolute irf weight
[c,i]=max(abs(IR));
lag = i-1;  % because i==1 corresponds to zero lag
if abs(c) >= abs(CL);
   s=1;
else
   s=0;
end

c = IR(i);
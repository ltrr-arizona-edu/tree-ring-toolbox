function h = highrow(L)
% highrow:  highest row with nonzero entry in logical matrix
% CALL: h = highrow(L);
%
% Meko -24-97
%
%***********************  IN *********************************
%
% L (mL x nL)L logical matrix
%
%***************************  OUT ***************************
%
% h (1 x nL)i   highest row with nonzero entry in each col;  NaN if 
%   no non-zero elements in col
%
%******************* NOTES ************************************
%
% Written as utility for hydacc1.m to mark latest month of year
% with temperature below freezing


[mL,nL]=size(L);
dumcv = (1:mL)';
A = repmat(dumcv,1,nL); % matrix, same size as L;

B = A .* L; % 

[x,h]=max(B);

b = ~any(L);
if any(b);
   h(b)=NaN;
end



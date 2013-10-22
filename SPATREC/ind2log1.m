function L = ind2log1(I,n)
% ind2log: index pointer (NaN filled) to a logical pointer
% CALL: L=ind2log(I,n)
%
%****************** IN ARGS ****************
%
% I (index array)
% n (desired col size of L;  row size of L will be same as I)
%
%***************** OUT ARGS ********************************
%
% L (mL x nL)L logical pointer matrix corresponding to index pointer I
%
%******************** NOTES********************************
%
% Example: I=[1 2 5 NaN; 3 4 NaN NaN] 
%          n specified as 7
%    gives L = [1 1 0 0 1 0 0; 0 0 1 1 0 0 0]
%
%   If I is empty, L is empty
%   If n is smaller than the largest value in v, the col size of L
%     becomes equal to the largest value in v.  For example, if n was
%     specified as 4 in the above example, we see that the largest element
%     of v is 5, and L becomes
%     L=[ 1 1 0 0 1; 0 0 1 1 0];

switch nargin
case 1
   n=[];
case 2
otherwise
   error('Number of in args must be 1 or 2');
end

   

[mI,nI]=size(I);
L =zeros(mI,n);

for i = 1:mI;
	v = I(i,:);
	L1 = ~isnan(v);
	n1 = sum(L1);
	%vones = ones(1,n1);
	L(i,v(L1))=1;
end 
L=logical(L);


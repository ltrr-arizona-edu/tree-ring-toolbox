function L = ind2log(I,n)
%
% Index array (NaN-filled) to a logical array
%
%****************** IN ARGS ****************
%
% I (index arrary)
% n (desired col size of L;  row size of L will be same as I)


[mI,nI]=size(I);
L =zeros(mI,n);

for i = 1:mI;
	v = I(i,:);
	L1 = ~isnan(v);
	n1 = sum(L1);
	vones = ones(1,n1);
	L(i,v(L1))=vones;
end 



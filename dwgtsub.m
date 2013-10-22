%B=rand(3,7)
%BB=B;
[mB,nB]=size(B);


[Y,I]=sort(B')

for i=1:mB
	B(i,I(1:nB-nkeep,i)) = zeros(1,nB-nkeep);
end




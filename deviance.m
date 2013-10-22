function d=deviance(n)
% Similarity measure for a terminal node of a tree-based model
% D Meko 11-18-96
% See Venable and Ripley (1994), p. 333
nsum=sum(n);
p=n/nsum; % probability of each class

fsum=0;
for i = 1:length(n);
	if n(i)~=0
		f = n(i)*log(p(i));
	else
		f=0;
	end
	fsum=fsum+f;
end
d = -2*fsum;

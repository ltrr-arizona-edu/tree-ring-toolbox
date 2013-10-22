function nn=struc3(na,nb,nc,maxnp)

% Build 3-col number-of-parameters matrix for input to, say, oe1.m
%
% D. Meko 12-21-93
%
%
%****************  INPUT ARGS  ***********************
%
% na - rv of number of first parameters. Like (1:4)
% nb - rv of no. of second parameters
% nc - rv of no. of third parameters
% maxnp - maximum allowable no. of total params in a model
%
%
%****************  OUTPUT ARS *******************************
%
% nn (? x 3) the num-of-params matrix
%
%
%************ NOTES *****************************************
%
% Used to prepare parameter number matrices for calls to oe1.m,
% arx1.m.  
% First 3 input args corresp to na, nb, nk for arx models,
%   to nb,nf,nk for oe models.


r=1;
for i=na;
	for j=nb;
		for k=nc;
			nn(r,:)=[i,j,k];
			r=r+1;
		
		end
	end
end
s=(sum(nn'))';
L=s<=maxnp;
nn=nn(L,:);
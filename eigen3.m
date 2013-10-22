function commun= eigen3(S,n)

% Compute communality from structure matrix, assuming n components used

% Feb 29, 1992, by D. Meko

%**********   INPUT ARGS
%
% S  (p x p)  structure matrix;  from running eigen1.m
% n  (1 x 1)  number of eigenvectors selected


%*******   OUTPUT ARGS
%
% commun (p x 1)  communality for each of p variables.
%	Equivalent to total pct variance expld for the variable by the
%	n selected components


[p,p]=size(S);

S(:,n+1:p)=[];  %  get rid of all factors beyond the n th from the
	% structure matrix

commun=diag(S * S');

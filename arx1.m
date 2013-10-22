function [V,iv,mod,th]= arx1(z,nn,n)

% Best ARX models for a time series 
%
% D. Meko 12-30-93
%
%
%**************** INPUT ***********************************
%
% z = [y u] (mz x2) time series, each row an observation (year)
%	col 1 is the output series, col 2 the input series; 
%	Usually these will have been detrended before being passed.
% nn (? x 3) orders na,nb,nk for each model (row)
% n (1 x 1) how many 'best" models to return in V
%
%************  OUTPUT ARGS **********************
%
% V (n x 2) Loss function (col 1) and FPE (col 2) of the
%	n models with lowest FPE; ranked #1 first
% iv (n x 1) row index to nn pointing to the n best models
% mod (n x 5) string vector of names of n best models


nnsize=size(nn);
for i=1:nnsize(1);
	th=arx(z,nn(i,:));
	VVV(i)=th(1,1);
	FPE(i)=th(2,1);
	mod(i,:)=[int2str(nn(i,1)),',',int2str(nn(i,2)),...
		',',int2str(nn(i,3))];
end


% Find the n best (lowest FPE models)
[fpe,ifpe]=sort(FPE);
fpe=fpe(1:n);
ifpe=ifpe(1:n);
VV=VVV(ifpe);
mod=mod(ifpe,:);

V=[VV;fpe]; % loss functions and final prediction error
iv=ifpe;


% Re-fit best model
nnn=nn(iv(1),:); % model order
th=arx(z,nnn);

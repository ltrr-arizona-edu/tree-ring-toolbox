function avcorr = rtree4(r)
%
% DESCRIPTION :  avcorr = rtree4(r)
% This function finds the average correlations between trees.
%
% INPUTS  :  r (? x 4) -  Correlation matrix. 
%		Column 1 - Correlation coeffs.
%		Column 2 - Sample size.
%		Coumns 3 & 4 - Sequence numbers of Tree 
%			ID's used to find the corr. coef.
%
% OUTPUTS  :  avcorr (1 x 2) - Average correlations and sample
%		size. 1st col. - Avg. corr. between trees.
%			2nd col. - Sample sizes.
%________________________________________________________________

% Initialize avcorr
avcorr=zeros(1,2);

L1 = ~isnan(r(:,1));
rx = r(L1,:);


% Average correlation between trees
avcorr(1,1)=mean(rx(:,1));
avcorr(1,2)=length(rx(:,1));

% End of file;

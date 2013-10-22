function avcorr = rtree2(r)

% DESCRIPTION :  avcorr = rtree2(r,cmask)
% This function finds the average correlations between trees
% and average correlations within trees.
%
% INPUTS  :  r (? x 5) -  Autocorrelation matrix. 
%		Column 1 - Correlation coeffs.
%		Column 2 - Sample size.
%		Column 3 - 1 if same tree, 2 if not.
%		Coumns 4 & 5 - Sequence numbers of Core 
%			ID's used to find the corr. coef.
%	   
%
% OUTPUTS  :  avcorr (2 x 2) - Average correlations and sample
%		size. 1st row - Avg. corr. between trees and 
%			within trees.
%			2nd row - Sample sizes.
%________________________________________________________________

% Initialize avcorr
avcorr=zeros(2,2);

% Average correlation between trees
r1=r(r(:,3)==2,1); 
if isempty(r1);
   avcorr(1,1)=NaN;
   avcorr(2,1)=NaN;
else;
   avcorr(1,1)=mean(r1);
   avcorr(2,1)=length(r1);
end;

% Average correlation within trees
r2=r(r(:,3)==1,1);
if isempty(r2);
   avcorr(1,2)=NaN;
   avcorr(2,2)=NaN;
else;
   avcorr(1,2)=mean(r2);
   avcorr(2,2)=length(r2);
end;

% End of file;

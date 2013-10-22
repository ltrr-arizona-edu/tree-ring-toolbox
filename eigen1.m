function [R,V,L,S,F,B]=eigen1(x,k);
% eigen1: eigenvectors, or principal components, of a data matrix
% [R,V,L,S,F,B]=eigen1(x,k);
% Last revised 2-29-92
%
% Eigenvector analysis of a data matrix. A common tree-ring application is 
% to orthogonalize a time series matrix of tree-ring indices in order to 
% reduce the dimensions of the tree-ring data before using them in a 
% climate reconstruction model
%
%*** INPUT ARGUMENTS 
%
% x (N x p)r  data matrix, p variables, N years
% k (1 x 1)i  options
%   k==1 pca on correlation matrix
%   k==2 pca on covariance matrix
%   k==3 pca on crossproducts matrix
%
%*** OUTPUT 
%
% R (p x p)r correlation (k=1), covariance (k=2) or cross-products 
%		(k=3) matrix of data x 
% V (p x p)r eigenvectors;  each column an eigenvector, ordered as in L
% L (p x p)r diag, eigenvalues, arranged from largest to smallest
% S (p x p)r factor structure.  S(j,k) gives correlation of j th time series with
%	k th PC-score
% F (N x p)r PC scores (see notes)
%
%*** REFERENCES
%
% Mardia, K., Kent, J., and Bibby, J., 1979, Multivariate Analysis: 
% Academic Press
%
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- statistics
%
%*** NOTES
%
% The eigenvectors, or PCs, V, are scaled eigenvectors -- scaled such that the
% sum of squares of each eigenvector (column of V) is 1.0.  The variable B
% stores the corresponding unscaled eigenvectors.
%
% Scores are computed by original data postmultiplied by the eigenvectors V(p x p).
% If k==1, the "data" are the z-scores of original time series;  if k==2, the "data" 
% is departure of the original series from their means;  if k==3, the data are the 
% untainted original data:  F = data * V
%
% Factor-score time series in cols of F do not generally have a std deviation of
% 1.0 (unless the the time series are converted to z-scores before calling eigen1.m).
% Fritts (1976) plots eigenvector 'amplitudes', which are equivalent to the PC scores.
%
% Some references refer to the scores F as having zero mean and unit standard deviation.
% We could get this version of the scores by computing scores as 
% F = data * B, where B is the (p x p) factor-structure matrix. The variable Fb is 
% stored in the program in case you want this "standardized" version of the scores.


%****************  Size arrays, compute means, z-scores
[N,p]=size(x);  % rows are years, cols are chrons
m=mean(x);  % rv of column means
x1=zeros(N,p);
z=zscore(x);  % compute z-scores of x
m1=m(ones(N,1),:);  % Dupe col means into N rows


%**********  Make cross-products, covariance or correlation matrix

if k==3;  % Use crossproducts matrix
	R = (1/N) * ( x' * x);  % Biased
elseif k==2; % Use covariance matrix (dispersion matrix)
	x1 = (x - m1);  %  subtract col means from x
	R = (1/N) *  (x1' *  x1); % Biased
elseif k==1;  % Use correlation matrix
	R=  (1/N )  *  (z' * z);
end

R = R * (N/(N-1));  % Convert to unbiased  version of R
%  If enabled, computed F or test data differs from C & L, p. 113.

%*************  COMPUTE EIGENVECTORS
% Re-order so that eigenvalues from largest to smallest

[V,L]=eig(R);   % eig is a matlab function to get eigenvectors V and
%		eigenvalues L of a square matrix R.

[yy,ii]=sort(diag(L));  % sort eigenvalues, smallest to largest
jj=ii(length(ii):-1:1);  % reverse order of eigvalues and vectors
V=V(:,jj);
L=diag(yy(length(yy):-1:1),0);


%******** COMPUTE MATRICES OF FACTOR-SCORE COEFFICIENTS; COMPUTE
%		FACTOR STRUCTURE

B = V / sqrt(L);  % Factor-score coefs;  factors are cols, sites are rows
	%  Division by sqrt (L) so that factor scores have mean zero and unit
	%  standard deviation.

S=  V * sqrt(L);  % Factor structure

if k==1; % PCA was on correl mtx
   Fb= z * B;  % non-standardized factor scores, version I do not use
   F = z * V; % factor score, version I use
elseif k==2; % pca was on covariance mtx
   Fb=x1 * B;
   F = x1 * V;
elseif k==3; %cross prod mtx
   Fb = x * B;
   F = x * V;
end


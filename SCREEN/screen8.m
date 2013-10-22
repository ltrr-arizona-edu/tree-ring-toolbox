function screen8
% screen9:  rotated pca exploratory analysis of screened tsm
% CALL:  screen8
% 
% Meko 8-19-97
%


% Get the time series matrix.  Assumes:
%   1. matrix name X
%   2. yr in col 1
%   3. 8-char series ids in N
[file1,path1]=uigetfile('tsm*.mat','Input file with time series matrix');
pf1=[path1 file1];
eval(['load ' pf1]);
if ~(exist('X')==1);
   error(['X missing from ' pf1]);
end
if ~(exist('N')==1);
   error(['N missing from ' pf1]);
end
if ~ischar(N);
   error('N must be char mtx');
end

% Pull year col off X
yr = X(:,1);
X(:,1)=[];

% Size X
[nyrs1,nvar1]=size(X); % number of years, number of variables

%Check rank of X to guard against duplicate columns
if rank(X)~=nvar1;
   error('Rank of X is less than number of variables cols in X');
end


% Call a NAG function for PCA, using correlation matrix
[E,P,V,S,IFAIL]= g03aaf(X,'C');
% E  (ith col):
%   1 eigenvalue assoc with ith PC
%   2 propo variation explained by ith PC
%   3 cum prop variation explained by first i PCs
%   4 Chi2 stat for ith PC
%   5 deg freedom for the chisq
%   6 (if method not 'C') sig level for the chi2 stat
% P the PC loadings; col j contains loadings for jth PC
% V the PC scores


%********************

% Stair plot of eigenvalues, with horiz line at 1.0
figure(1);
plot(E(:,1));
title('Eigenvalues')
line([1 nvar1],[1 1]);

% Scree plot, with horiz lines at .9 an .95
figure(2);
stairs(E(:,3));
set(gca,'Ylim',[0 1]);
title('Cumulative Variation Explained');
xlabel('PC #');
line([1 nvar1],[.90 .90]);
line([1 nvar1],[.95 .95]);

% Set threshold for number of PCs to retain for rotation
prompt={'Enter the number to retain:'};
def={sum(E(:,1)>1.0)};
title=['Retain this many of ' int2str(nvar1) ' PCs for Rotation'];
LineNo=1;
answer=inputdlg(prompt,title,LineNo,def);
nkeep = str2num(answer{1});

% Get proportion of variation explained by retained PCs
pkeep = E(nkeep,3);
 


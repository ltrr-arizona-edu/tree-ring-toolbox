function colbar01(colmap,ll,ur,txtlab,fontsz)
% colbar01:  colorbar of patches with labels at right
% CALL: colbar01(colmap,ur,txtlab,fontsz);
%
% Why?  First needed to label color map from resp02.m showing sign and 
% significance of largest impulse response weight of tree rings to 
% climate
%
% Meko 10-31-97
%
%*******************  IN ****************************
%
% ll (n1 x 2)r x,y coordinates of lower left part of each patch
% ur (n1 x 2)r x,y coords of upper right of each patch
% txtlab(n1 x ?)s label for each patch
% fontsz(1 x 1)r   font size for label 
% colmap(? x 3)i  triplets of colors for the n1 points defining a color map
%
%******************** OUT ****************
%
% No out args
%
% Places a color legend on a current map.  For example, might
% want 4 rectangles of different colors, indicating whether
% weight strong and positive, weak and positive, strong and negative,
% or weak and negative
%
%*************  END OF OPENING COMMENTS


%***********  CHECK INPUT
[mC,nC]=size(colmap);
if nC ~=3;
   error('colmap must be 3-col');
end


colormap(colmap);

[m1,n1]=size(ll);
npatches=m1;
[m2,n2]=size(txtlab);
[m3,n3]=size(ur);
if n3~=2;
	error('ur must be 2-col matrix');
end


Ltemp=[m2==m1  m3==m1 mC==m1];
if ~all(Ltemp);
	error('colmap,ll,ur,txtlab, must be same row size');
end


% Compute four vertices for each of n1 patches
ul = [ll(:,1) ur(:,2)]; % upper left
lr = [ur(:,1) ll(:,2)]; % lower right

% Build X matrix
X = [(ll(:,1))'; (ul(:,1))'; (ur(:,1))';  (lr(:,1))'];

% Build Y matrix
Y = [(ll(:,2))'; (ul(:,2))'; (ur(:,2))';  (lr(:,2))'];


% Compute plotting point for left edge of labels
incx = 0.10* (ur(1,1)-ll(1,1));
mr = [incx+lr(:,1)      (lr(:,2) + ur(:,2))/2];



% ************  Plot the patches


patch(X,Y,[1 2]);



% ************  Add labels
for n = 1:npatches;
	xtxt = mr(n,1);
	ytxt = mr(n,2);
	text(xtxt,ytxt,txtlab(n,:),'FontSize',fontsz);
end












%function dat2mat1(file1,nch,wch,sfx)
% Convert ".sfx" files to mat files with same prefix
% D. Meko 10-24-94
%
%
%
%****************** IN ARGS *********************
%
% file1 -- name of file (in current directory)  holding string
%   matrix with individual station file names
% nch -- col size of the string matrix
% wch (1 x 2) -- which cols hold the  individual filename prefixes
% sfx -- string suffix (without the dot) assumed for all the station files
%
%


B = txtin1(file1,nwh);
[mB,nB]=size(B);
for i=1:mB;
	str = B(i,wch(1):wch(2));
	eval(['load ',str,'.',sfx]);
	eval(['X = ',str,';']);
	eval(['save ',str,' X']);
end






% A (mA x nch) 

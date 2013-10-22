function T=txtin1(fn,nc)
% Read in a matix of text strings
% D. Meko;  10-12-94
%
%************************************   INPUT
%
% fn (1 x ?)s   file name, including path if needed
% nc (1 x 1)i   number of characters in each row of the string matrix
%
%
%************** OUT ARGS
%
% T (? x nc)s   text matrix
%
%
%*************   NOTES **********************************
%
% Use the character "@" in place of blank spaces in the ascii file fn
%
%*************************************************************

if nargin ~=2
	error('Number of args not 2')
end


[fid,message] = fopen (fn,'rt');
if fid == -1;
	disp(['Message from fopen: ',message])
	disp(['fid = ',int2str(fid),':  no luck opening ',fn])
	error('Check that the file exists and is in path')
end
[f,count] = fscanf(fid,'%s');
disp(['Number of ',int2str(nc),'-char elements read = ',int2str(count)])
f = f';  % convert to a cv
nlen = length(f);

% Check that length of f consistent with specified number of columns
n1 = rem(nlen,nc);
if all(n1)~=0;  % 
	error(['String length ',int2str(nlen),...
	' not evenly divisible by # of cols ',int2str(nc)]);
end
nr = nlen/nc;  % Number of rows in the string matrix

f = reshape(f,nc,nr);
f = f';

% Replace the @ with blank
L = f == '@';
nL = sum(sum(L));
blank = ' ';
f(L) = blank(ones(nL,1),:);
T=f;

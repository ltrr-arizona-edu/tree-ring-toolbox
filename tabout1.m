function tabout1(D,fmt,fnm)
% tabout1: matrix to a tab-delimited output file, without column headings
% CALL: tabout1(D,fmt,fnm);
%
% D Meko 12-19-95
%
%************** IN ARGS *****************
%
% D (mD x nD) -- data matrix
% fmt (s) -- format for printing  
% fnm (s) -- filename for output
%
%*************** NOTES ****************
%
% If only two input args, you will be prompted for the filename
% of the output file. 
%
% Example of fmt: '%5.2f\t%5.2f\t%5.3f\n'  -- for ? x 3 matrix
% \t puts tabs between cols, as required if you want to open the
% file in MS Word to convert to a table.  \n puts each row on a
% new line
%
% Input format takes care of rounding and truncating
%
% User can open the ascii output file in WordPerfect or MS word.  Can also insert
% the file into document



% Size data matrix
[mD,nD]=size(D);
mn = mD*nD;  % total number of elements in D

% Reshape the matrix into a col vector by transposing and then
% stacking colums
D = D';
D=reshape(D,mn,1);

if nargin==2; % user wants to be prompted for output file name
	fnm=input('Output Filename?','s')
end

% Next line for debugging only -- writes to screen instead of
% file
%fprintf(1,fmt,D)


% Create the output file with write as text permission, or
% open the file if it exists -- it will be overwritten
fid = fopen(fnm,'wt')

% write the data into the file
fprintf(fid,fmt,D);

% Close the file
fclose(fid)


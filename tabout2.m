function tabout2(D,H,fmtd,fmth,fnm)
% tabout2: mtx to tab-delimited file; with col headers
% CALL: tabout2(D,H,fmtd,fmth,fnm);
%
% D Meko 12-20-95
%
%************** IN ARGS *****************
%
% D (mD x nD) -- data matrix
% H (mH x nH)s each row is a header for col of D
% fmtd (s) -- format for printing  D
% fmth (s) -- format for printing one header
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
% Example of fmtd: '%3s' -- for 3-char header
%   tabout2.m will handle putting on the \t and \n as needed
%
% Input format takes care of rounding and
% truncating.


% Size data matrix
[mD,nD]=size(D);
mn = mD*nD;  % total number of elements in D

% Size header matrix
[mH,nH]=size(H);
if mH~=nD, error('row size of H must equal col size of D'), end;

% Reshape data matrix into a col vector by transposing and then
% stacking colums
D = D';
D=reshape(D,mn,1);

if nargin==4; % user wants to be prompted for output file name
	fnm=input('Output Filename?','s')
elseif nargin==5; 
	% user had put file name in arg list'
else
	error('number of input args should be 4 or 5')
end

% Next line for debugging only -- writes to screen instead of
% file
%fprintf(1,fmt,D)


% Create the output file with write as text permission, or
% open the file if it exists -- it will be overwritten
fid = fopen(fnm,'wt')

% write the headers for columns
for i=1:mH;
	h = (H(i,:))';
	if i==mH,
		fmthh=[fmth '\n'];
	else
		fmthh=[fmth '\t'];
	end
	fprintf(fid,fmthh,h);
end


% write the data into the file
fprintf(fid,fmtd,D);

% Close the file
fclose(fid)


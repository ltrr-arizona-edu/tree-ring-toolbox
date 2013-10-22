function [H,X]=loadhead(flnm,ncols);
% loadhead: read in time series matrix with column headings using fscanf; BMDP 
% CALL: [H,X]=loadhead(flnm,ncols);
%
% Meko 6-29-97
%
%
%***************  IN ***************************
%
% flnm (1 x ?)s filename, possibly including path of file with input data
% ncols (1 x 1)i number of columns, in the input matrix, including years col
%
%**************** OUT *******************************
%
% H (nsers x  nchars)s   matrix of headings for the cols 
%			nsers is one less than ncols
%			nchars is the maximum length of any heading input
%
%
%******************** NOTE ********************************
%
% Some heading for the years column is assumed in the input file
% Column headings are delineated by spaces

hmsg=msgbox('NOTE: THIS FUNCTION USES FSCANF, WHICH WILL NOT WORK WITH NaNs');


% Open input file
fid1=fopen(flnm,'r');




%******************* BUILD STRING MATRIX OF COLUMN HEADINGS

% Get the header line
h = fgetl(fid1);
hlen = length(h);

h = strjust(h); % right justify the header

% Remove any leading blanks
h=fliplr(deblank(fliplr(h)));

% h should now begin with the first char of the year header and end with
% the last char of the last col header.

% Add a single space at left of header and right.  Then any change from space to 
% non space will mark start of a col header; any change from 
% nonspace to space will mark end of header
h = [ ' ' h ' '];

hs= isspace(h);
dh = diff(hs);
igo = find(dh== -1) + 1; 
isp = find(dh==1) + 1;

% Check that number of identified headers equals number of cols
if (length(igo) ~= length(isp)) | length(igo)~=ncols;
   error('Number of identified column headings does not match number of cols');
end

% Find max  number of chars in any col heading, and initialize storage matrix for names
max1 = max(isp-igo+1);
H = blanks(max1);
H=repmat(H,ncols,1);

% Store headings, left justified
ilen=isp-igo+1;
for n = 1:ncols;
   i1=igo(n); i2=isp(n);
   H(n,1:ilen(n)) = h(i1:i2);
end



%**********************  GET THE DATA FOLLOWING THE HEADINGS

% Put data in col vector
x=fscanf(fid1,'%g');

% Check that length x divisible by number of cols
if rem(length(x),ncols)~=0;
   error('Length of x does not equal ncols');
end

% Compute number of years of data
nyrs = length(x)/ncols;

% Re-form the matrix
X=(reshape(x,ncols,nyrs))';
fclose(fid1);

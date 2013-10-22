function datfrm1(X,flnm,fmt,hdcol,hdrow)
% datfrm1:  convert matrix into an SPLUS dataframe, with optional col names and row names
% CALL: datfrm1(X,flnm,fmt,hdcol,hdrow);
%
% Meko 5-26-98
%
%************ IN 
%
% X (mX x nX)r  all-numeric matrix to be converted
% flnm (1 x ?)c  path\filename for output dataframe
% fmt (1 x ?)c  format for writing  a row of data to dataframe
% hdcol {1 x nX}c column headings----OPTIONAL
% hdrow {mX x 1)c row headings--------OPTIONAL, and only possible if also have hdcol
%
%*********** OUT
%
% No output args
%
% Output file with path/name   flnm holds the ascii dataframe 
%
%************** NOTES
%
% User of SPLUS can import the file in flnm as a dataframe. In SPLUS, settings
% for import will depend on whether
%  no col headings or row headins
%  col headings but no row names
%  col headings and row names
%
% You can have no col names or row names, col names but no row names, or col names and
% row names.  You cannot have row names and no col names.
%
% In SPLUS, will use file/import data/from file/options...
% Set namerow to 1 if have col headings in hdcol;  set namecol to 1 if have row
% names in hdrow.  Otherwise, leave namerow and namecol blank in the input dialog



%---------- CHECK INPUT

[mX,nX]=size(X);

if ~ischar(flnm) | size(flnm,1)~=1; 
   error('Filename must be char row vector');
end

if ~ischar(fmt) | size(fmt,1)~=1; 
   error('fmt must be char row vector');
end

if nargin>3;
   if ~iscell(hdcol);
      error('hdcol must be cell');
   end
   if size(hdcol,1)~=1;
      error('hdcol must have 1 row');
   end
   if size(hdcol,2)~=nX;
      error('col size of X must match col size of hdcol');
   end
end
if nargin>4;
   if ~iscell(hdrow);
      error('hdrow must be cell');
   end
   if size(hdcol,1)~=1;
      error('hdrow must have 1 row');
   end
   if size(hdrow,1)~=mX;
      error('col size of X must match col size of hdrow');
   end
end
     

%---------------------- Handle column headers

if nargin>3;
   
   % Make sure no internal spaces in column headings
   for n = 1:nX;
      cc=  char(hdcol(1,n));
      if any(isspace(cc));
         error('a col heading has an internal space');
      end
   end
   
   hdcol2 = char(hdcol); % convert cell to character matrix, one label per line
   [m1,n1]=size(hdcol2);
   hdcol2 = [hdcol2 repmat(blanks(1),m1,1)];  % add blank to right hand side
   hdcol2= hdcol2'; % transpose, because will write down cols
   hdstrcol = sprintf('%s',hdcol2);
   if nargin==5;
      hdstrcol = ['rownames ' hdstrcol];
   end
   
   
end


%---------------------- Handle row headers

if nargin==5;
   
   % Make sure no internal spaces in column headings
   for n = 1:mX;
      cc=  char(hdrow(n,1));
      if any(isspace(cc));
         error('a row heading has an internal space');
      end
   end
   
   hdrow2 = char(hdrow); % convert cell to character matrix, one label per line
   [m1,n1]=size(hdrow2);
   hdrow2 = [hdrow2 repmat(blanks(1),m1,1)];  % add blank to right hand side
end


%-------------------- Open out file
fid1 = fopen(flnm,'w');



%--------------------- Write dataframe

if nargin==3; % No col or row headings
   fprintf(fid1,fmt,X');


elseif nargin==4; % Col headings, but no row headings
   fprintf(fid1,'%s\n',hdstrcol);
   
   %------ write numeric data matrix
   fprintf(fid1,fmt,X');
elseif nargin==5;  % col headings and row headings
   fprintf(fid1,'%s\n',hdstrcol); % write col headings
   for n = 1:mX;
      thisrow = hdrow2(n,:);
      x = X(n,:);
      fprintf(fid1,'%s',thisrow);
      fprintf(fid1,fmt,x');
   end
else;
   error('valid argin is 3,4, or 5');
end

fclose(fid1);


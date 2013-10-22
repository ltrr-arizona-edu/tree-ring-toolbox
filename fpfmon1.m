function fpfmon1(Y,fn,tit,digs,dtype)
% Print monthly ppt matrix for appendix of report
% CALL: fpfmon1(Y,fn,tit,digs,dtype)
%
% D Meko 10-27-94
%   Rev: 9-25-96
%
%********************  IN ARGS  ***************************
%
% Y (mY x 13)r   year in col 1, 12 ppt values in inches for jan-dec
% fn (?)s  name for output ascii file (e.g., 'TUCSON.TXT')
% tit (1 x 64)s  title 
% digs (1 x 1)i scaling digit, typically 2 for ppt data
%    Data is multiplied by 1E"digs", then rounded to nearest integer;
%    Digs=2 gives 1E2=100 as scaling, which produces hundredths of 
%	  inches if Y is input as inches
% dtype
%	1=pcp
%	2=tmp
%***************************  OUT ARGS -- NONE
%
%********************  NOTES  **********************
% 
% Scales and rounds so that printed data to hundredths of inch, 
%  no decimal point.
% Year and month headings hard-coded
% Numeric field for pcp is 5, allowing for max of 99.99 inches (9999) with
%  one blank spacing between months


% Check input mtx
[mY,nY]=size(Y);
if(nY~=13)
	error('ppt matrix should be 13 cols')
end

% Scale and round so that ascii file is in hundredths of inches;
% Transpose because fprintf will grab data down cols
tenfact=eval(['1e',num2str(digs)]);  % scaling factor; 
% dig=2 gives 100, the usual for TRL ppt data
Z = round([Y(:,1) (Y(:,2:13)*tenfact)]);

% Build format for a line of printout
if dtype==1;
	% Build format for a line of printout if data is pcp
	fmt1='%4.0f%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f';
	fmt2='%5.0f%5.0f%5.0f%5.0f%5.0f%5.0f\n';
elseif dtype==2; % tmp data
	% Build format for a line of printout if data is tmp
	fmt1='%4.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f';
	fmt2='%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f\n';
else
	error('dtype must be 1')
end
	fmt=[fmt1 fmt2];


% open file output file for writing
fid=fopen(fn,'w');

% Build title
blnk=blanks(64);  % make 1 x 64 string of blanks
lt=length(tit); % Number of characters in string argument for title
if lt > 64
	error('Make title 64 or fewer characters')
elseif lt==64
	tit1=tit;  % use title as is
else;  % center title
	ndiff=64-lt;
	nn=fix(ndiff/2);
	remmy = rem(ndiff,2);
	if remmy==0
		tit1 = [blanks(nn)  tit blanks(nn)];
	else
		tit1=  [blanks(nn+remmy) tit blanks(nn)];
	end
end
%end

% Print title followed by blank line 
fprintf(fid,'%s',tit1);
fprintf(fid,'\n\n');


% Print header line; after checking that 
hdr = 'Year   J    F    M    A    M    J    J    A    S    O    N    D';

%hdr1= 'Year     J     F     M     A     M     J     J';
%hdr2= '     A     S     O     N     D';
%hdr=[hdr1 hdr2];

fprintf(fid,'%s',hdr);
fprintf(fid,'\n\n');  % position at next line for data


for i=1:mY;    % loop over years
	fprintf(fid,fmt,Z(i,:)');
end


fclose all;


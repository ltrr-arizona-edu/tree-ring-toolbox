function rwwrite(X,person,when,fln)
% rwwrite(X,person,when,fln)
%
% Write ringwidths to a file with name "fln"
%
% D Meko 12-31-95
%
%****************** IN ARGS 
%
% X (mX x nX)r measurements in col 2, year in col 1
% person (1 x 3)s initials of measurer
% when (1 x 8)s day/month/year measurements recorded
% fln ()s name of output file containing ringwidths, 
%    ".rw" not included
%
%***************** OUT ARGS -- NONE
%
%***************** UW FUNCTIONS CALLED -- NONE
%
%***************** NOTES
%
% Needed as part of rwlook.m function group to generate '.rw' 
% file formatted output of edited ring-width data that could
% be read by convert.exe (for .rwl file from .rw files)
%
% X is generated in the calling program. person and when are
% lines 1 and 2 headers from .rw file.



% X must be 2-column matrix, with elements in first column increasing
% by increments of 1 down the column (year or nominal year), and
% elements of column 2 non-negative and less than 999.
[mX,nX]=size(X);
if nX~=2;
	error('X must be 2-col matrix')
end
d1=diff(X(:,1));
d2=d1==1;
if ~all(d2),
	error ('Column 1 of X must increment by 1')
end
L1=X(:,2)<0 | X(:,2)>=999;
if any(L1),
	error('Measurements in col 2 of X must be + and  <999')
end

% truncate or expand the measurer's name to 3 characters, if necessary
if length(person)>3; init=person(1:3); end
if length(person)<3;
	npad=3-length(person);
	init=[person blanks(npad)];
end
if length(person)==3,
	init=person;
end

w=when(1:8);

% All systems go -- open the file for writing
fid=fopen(fln,'wt'); 

% Write the measurer's initials, measurement date, and
% first year of measurements
fprintf(fid,'%3s\n',init)
fprintf(fid,'%8s\n',w)
fprintf(fid,'%5.0f\n',X(1,1))

% Put measurements in a column vector and write to file
x=X(:,2);
fprintf(fid,'%5.0f\n',x)

% Put 999 to mark end of measurements
iend=999;
fprintf(fid,'%5.0f',iend)

fclose(fid)

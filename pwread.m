function [X,fln]=rwread

% Reads a file of data in ".rw" format
% Puts year in col 1, ring width in col 2

% X (mX x 2) year in col 1, ring width tin col 2
% fln ( name of file containing ring width (rw file)
%     Be sure path includes directory containing file.

% D. Meko 10-27-93

%*******  USER-WRITTEN FUNCTIONS  -- NONE

global fln;  % must be able to find the fil in global space

% prompt for rw file name.  Then add suffix .rw
fln = input('Name of RW file (without suffix): ','s');
fln = [fln '.rw'];

% Open rw file for reading
fid = fopen(fln,'r');


% Read and echo first 2 lines of file.  Line 1 is the measurer;s
% initials.  Line 2 is the date the core was measured.  
disp(['Measured by  '  fgetl(fid)  ' on '  fgets(fid)]);

% Read the beginning year for measurements
yrgo = fscanf(fid,'%d',1);

% Read measurements into a cv;
x=fscanf(fid,'%d',inf);
len=length(x);  % How many years of ring width? Still icludes the
	% final 999 value
yr=(yrgo:yrgo+len-1)';  % cv of years

% Delete any 999 rows -- usually just the last row
L=x==999;
x(L,:)=[];
yr(L,:)=[];
disp([num2str(sum(L))  ' 999 values(s) were found in the rw file']);
disp('press RETURN to continue');
pause(2)

% Put years in col 1 , measurements in col 2 of X
X=[yr x];

fclose(fid)

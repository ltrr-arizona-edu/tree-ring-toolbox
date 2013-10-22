function [X,person,when,file1]=rwread(pre,txt1)
% Read in a ring-width series in .rw format for use by rwlook.m, for example
% CALL:  [X,person,when,file1]=rwread(pre,txt1);
%
% D. Meko 10-27-93, 9-20-96 to use uigetfile function 

% Reads a file of data in ".rw" format
% Puts year in col 1, ring width in col 2
%
%
%*************  IN
%
% pre -- cell array of .rw file prefixes
% txt1 -- string, either "first", or "second", indicating which  series
%   user is prompted to select
%
%************* OUT
%
% X (mX x 2) year in col 1, ring width tin col 2
% fln ( name of file containing ring width (rw file)
%     Be sure path includes directory containing file.


%*******  USER-WRITTEN FUNCTIONS  -- NONE

global fln;  % must be able to find the fil in global space

% Get a .rw file
%[file1,path1]=uigetfile('*.rw','Select ring-width file');
%pf1=[path1 file1];


kmen1=menu(['Choose ' txt1 ' .rw file'],pre);
file1 = char(pre(kmen1));
file1 = [file1 '.rw'];
pf1 = file1;

% Open rw file for reading
fid = fopen(pf1,'r');

person=fgetl(fid); % name of measurer
when=fgets(fid);

% Read and echo first 2 lines of file.  Line 1 is the measurer;s
% initials.  Line 2 is the date the core was measured.  
%disp(['Measured by  '  fgetl(fid)  ' on '  fgets(fid)]);
disp(['Measured by  '  person  ' on '  strtok(when)]);

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

n999 = sum(L); % number of 999 values found in files
if n999==1;
   disp([num2str(sum(L))  ' 999 value was found in the rw file, as expected']);
   pause(1);
else
   disp([num2str(sum(L))  ' 999 values were found in the rw file']);
   pause(3);
end

% Put years in col 1 , measurements in col 2 of X
X=[yr x];

fclose(fid);

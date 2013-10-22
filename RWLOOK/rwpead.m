function [X,guy,day,fln]=rwpead(k1)
% [X,guy,day,fln]=rwpead(k1)
% Reads a file of data in ".rw" format
% Puts year in col 1, ring width in col 2
% Differs from rwsead in that uses uigetfile
%
% D. Meko 6-7-96
%
%******************** IN ARGS ***********************
%
% k1 -- 1=eww,  2=lww, 3=rw
%
%
%******************* OUT ARGS ************************
%
%
% X (mX x 2) year in col 1, ring width tin col 2
% guy -- string, initials of measurer
% day -- string, day measured
% fln ( name of file containing ring width (rw file)
%     Be sure path includes directory containing file.

%*******  USER-WRITTEN FUNCTIONS  -- NONE


if k1==1; % earlywood file
	ftype='*.eww';
	txt='Earlywood filename: ';
elseif k1==2,
	ftype='*.lww';
	txt='Latewood filename: ';
elseif k1==3;
	ftype='*.rw';
	txt='Total width filename: ';
else
	error('allowable responses 1,2, or 3')
end

% Prompt for the input file
[fln,pth]=uigetfile(ftype,txt);
pf1=[pth fln];

% Open rw file for reading
fid = fopen(pf1,'r');

guy=fgetl(fid); % name of measurer
day=fgets(fid);

% Read and echo first 2 lines of file.  Line 1 is the measurer;s
% initials.  Line 2 is the date the core was measured.  
%disp(['Measured by  '  fgetl(fid)  ' on '  fgets(fid)]);
disp(['Measured by  '  guy  ' on '  day]);

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
if sum(L)==1 & L(length(L))==1,
	disp('As expected, series has one 999 value')
	disp('and that is the last value in the series')
	pause(1)
elseif sum(L)==0
	disp('No trailing 999 value in series -- might be problem')
	disp('press RETURN to continue');
	pause
elseif sum(L)==1 & L(length(L))~=1,
	error('Single 999 value in series was not last value')
else
	error('More than one 999 value read')
end


% Put years in col 1 , measurements in col 2 of X
X=[yr x];

fclose(fid);

function X=rwread1(k1,root,epath,wpath)
% X=rwread1(k1,root,epath,wpath)
%
% Reads .rw, .eww, or .lww files for use in function ellook.m
%
% D. Meko 10-14-96
%
%******************** IN ARGS ***********************
%
% k1 -- 1=eww,  2=lww, 3=rw
% root -- root of the filename. Like: 'dcy04'
% epath -- path to .eww and .lww files (early and late growth)
% wpath -- path to .rw files (total ringwidth)
%
%******************* OUT ARGS ************************
%
% X (mX x 2) year in col 1, ring width tin col 2
%
%*******  USER-WRITTEN FUNCTIONS  -- NONE


if k1==1; % earlywood file
	ftype='.eww';
	fpath=epath;
	txt='Earlywood filename: ';
	root=[root 'xe'];
elseif k1==2,
	ftype='.lww';
	fpath=epath;
	txt='Latewood filename: ';
	root=[root 'xl'];
elseif k1==3;
	ftype='.rw';
	fpath=wpath;
	txt='Total width filename: ';
else
	error('allowable responses 1,2, or 3')
end

% Build file name
pf1=[fpath root ftype];

% Open rw file for reading
fid = fopen(pf1,'r');

guy=fgetl(fid); % name of measurer
day=fgets(fid);

% Read and echo first 2 lines of file.  Line 1 is the measurer;s
% initials.  Line 2 is the date the core was measured.  
%disp(['Measured by  '  fgetl(fid)  ' on '  fgets(fid)]);
disp(['Root of filename:  '  root]);
disp(['File type = ' ftype]);
disp(['Measured by  '  guy  ' on '  day]);

% Read the beginning year for measurements
yrgo = fscanf(fid,'%d',1);

% Read measurements into a cv;
x=fscanf(fid,'%d',inf);
len=length(x);  % How many years of ring width? Still includes the
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
function clicon4
% pptcon4: convert garfin-style set of hcn ppt or tmp files to .mat files
%
% Meko 5-8-97
%
%
%**************  INPUT ****************************
%
%
% User is prompted for the following:
%
%  flnms*.txt -- list of files with the monthly ppt data.  First line 
%			give source directory.  Second line gives target directory for the
%			output .mat files.  Remaining lines list the file suffixes.  
%
%			Example:
%
%				d:\greggin\
%				d:\greggout\
%				tucs.dat
%				bisb.dat
%				doug.dat
%				...etc
%
%			The trailing sequential number is not used by the program. Just 
%			convenient.
%
%
% The files tucs.dat, bisb.dat, etc are ascii data files, 13 cols, with
% a year and 12 values of monthly ppt per line.  Properties of this file:
%
%  - file is all numeric
%  - may have entire years' gaps
%  - missing monthly values coded by a number, which user is prompted for in 
%		a dialog box
%  - data unit might be different than desired in output. For example,
%		input might be  hundredths of inches of ppt, and output might be mm.
%		Conversion specified in dialog boxes.
%  
%
%************************ OUT *****************************
%
% * Set of .mat files, one per station, with year and monthly data in X
%   Data has these properties:
%
%		-data units might have been converted, depending on user responses to
%		 prompts
%		-any missing gaps of years have been replaced with the year and 12 NaN
%		-missing monthly values are NaN
%
% * Ascii file of the list of filenames of .mat files generated.  Truncated
%   prefix to 4 chars.  No '.' or suffix. Sequential number to right
%	  of filename.  For example:
%
%		tucs 1
%		bisb 2
% 	 etc


%-------------- GET LIST OF FILES

[file1,path1]=uigetfile('flnms*.txt','Input list of filenames');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

% Get source and target paths from lines 1 and 2
c=fgetl(fid1);
if c(2)~=':';
	fclose all
	error('Line 1 is not valid path');
else
	path2=strtok(c); % path to input monthly data files
end

c=fgetl(fid1);
if c(2)~=':';
	fclose all
	error('Line 2 is not valid path');
else
	path3=strtok(c); % path to output .mat files
end


%-----------------  INTIAL READ TO COUNT NUMBER OF MONTHLY DATA FILES

nfiles=0;
k1=1;
while k1;
	c=fgetl(fid1);
	if ~feof(fid1)
		if any(c=='.'); % probably a valid file name
			nfiles =nfiles+1;
		else; % might be and eof on an empty line
			k1=0;
		end
	else;  % reached an eof
		k1=0;
	end
end
clc
disp(['Number of data filenames counted = ' int2str(nfiles)]);


% rewind and position at first filename
frewind(fid1);
c=fgetl(fid1);
c=fgetl(fid1);


%------- Initialize

% String matrix to hold names
S1=blanks(12);


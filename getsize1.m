function [n,YRS]=getsize1(pf1)
%
% For a group of .mat pcp or tmp files, and a file with the names of
% those files,  get the number of files, and the start and end year
% of each
%
% Meko, 8-3-96
%
%
%*************** IN ARGS
%
% pf1 (?)s  path/filename of file with names of the .mat files
%		pf1 would typically be set in the calling function using 
%		uigetfile
%
%*************** OUT ARGS
%
% n (1 x 1)i  number of .mat files in pf1
% YRS (? x 2)i start end years of the n .mat files 
%
%
%****************** NOTES
%
% First used in SRP study to get info for setting size of matrix
% for seasonalized pcp data by seaspt4.m


% Allocate
a=NaN;
YRS = a(ones(1000,1),ones(2,1)); % hard code 1000 as max #
	% of files will ever work with
%warndlg('YRS is hard coded to handle up to 1000 .mat file',...
%	'Notice');
%warddlg('Path to .mat data files is hard coded-- see path2',...
%	'Notice');


% Open the file of filenames
fid1=fopen(pf1,'r');


n=0; % counter for number of files
k1=1;
while k1;
	c=fgetl(fid1);
	if feof(fid1);
		k1=0;
	else
		n=n+1; % increment file counter

		% Get the non-white-space part of file name
		c=strtok(c);
		c=c(~isspace(c));

		% Build the path/filename
		path2='d:\srp\ppt51\'; % hard coded path for data files
		pf2=[path2 c];
		disp(['Loading: ' pf2])

		% Load the .mat file and get its start and end year; note
		% that I usually store the data matrix as Z, but sometimes
		% as Y1, so need to check for both. Will store data in X
		eval(['load ' pf2]);
		if exist('Z') & ~exist('Y1')
			X=Z;
		elseif exist('Y1') & ~exist('Z')
			X=Y1;
		elseif exist('Y1') & exist('Z')
			error([' Both Z and Y1 found in file # '  int2str(n)])
		else
			error([' Neither Z nor Y1 found in file# ' int2str(n)])
		end

		
		% Store start, end year for this series in YRS
		YRS(n,:)=[min(X(:,1))  max(X(:,1))];
	end
end

	% Truncate YRS
	YRS=YRS(1:n,:);



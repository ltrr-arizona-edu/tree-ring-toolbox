function [Y,yr]=seaspt4(months,dtype,O,nrmpd)
% [Y,yr]=seaspt4(months,dtype,O,nrmpd) 
% Seasonalize monthly pcp or tmp data stored in .mat files
% Store the seasonalized data in a tsm, each col for a station
% 
% Meko 8-3-96
%
%*************************** IN ARGS *********************
%
% months (1 x 2)i start, end months of "season"
%		1=Jan, 2=Feb,... 12=Dec
%		ex: [11 4]  is "nov-apr"
% dtype (1 x 1)i data-type
%		1=pcp
%		2=tmp
% O (1 x 2)i options
%	O(1): real units or %/departure units for output
%		1==real units (e.g., in. for pct, deg F for tmp)
%		2==convert pcp to % of normal, tmp to departure from norm
%	O(2): normalizing period method
%		1==use entire available data for each series to compute normal
%		2==specify uniform normalizing period for all series
% nrmpd (1 x 2)i normalization period (optional; need only if
%		O(2)==2
%
%
%*************************** IN  FILES **************************
%
% User is prompted for name of file containing filenames of the
% .mat files for all stations;  the path/filename is stored in pf1
%
%****************** OUT: CONTENTS OF .MAT FILE
%
% User is prompted for name of a ".mat" storage file that contains
% key output. The stored matrices are listed below:
%
% Y (mY x nY)r  the seasonalized data; mY years and nY stations
% yr (mY x 1)i year vector for Y
% YRS1 (nY x 2)i  start, end years of .mat files input with
%		monthly data
% YRS2 (nY x 2)i  start, end years of valid (non-NaN) seasonalized
%		data for each series;
% S (2 x nY)r  means and sample sizes for normalizing
%   Period for computation depends on O(2)
%		row 1: means 
%		row 2: sample size (number of years) for conputing normals
% months, dtype, O, nrmpd,pf1 :    see input
%*********************** NOTES *****************************
%
% First used in SRP 1996 study
% Seasonal tsm automatically col-sized to accomodate however many
%	stations are in the station-names file, and row-sized according
%	to earliest start year and latest end year of any monthly series
% Output can be in same units as monthly data, 
%	or in pct normal (for pcp) or departure from normal units.
% If in pct or departure units, normal period may be all available
%	valid seasonal values or the series, or a specified normalizing
%	period the same for all series.  If a specified normalizing 
%	period, however many years of data are available for that period
%	will be used in computing the normal, and the user will be
%	warned if any series does not have complete coverage for the 
%	normal period
% 


a=NaN;

% Set path/filename to get names of mat files with monthly data
[file1,path1]=uigetfile('*.txt','File with filenames of monthly data');
pf1=[path1 file1];

% Compute number of montly-data files and their start, end years
[n,YRS1]=getsize1(pf1);

% Allocate Y to store seasonalized data; initialize as NaN
% Vector yr will hold the years for Y
yr1=min(YRS1(:,1));
yr2=max(YRS1(:,2));
yr = (yr1:yr2)';
nyrs=length(yr);
Y = a(ones(nyrs,1),ones(n,1)); 

% Allocate storage for means and sample sizes
S=a(ones(2,1),ones(n,1));


% Open the file of filenames of .mat monthly dat files
% Already know that there are n of these files from
% call to getsize1.m
fid1=fopen(pf1,'r');

% Loop over files, getting their monthly data and seasonalizing
clc
disp('Starting Computation of Seasonalized Data')
for nn = 1:n;
	c=fgetl(fid1);
	c=strtok(c);
	c=c(~isspace(c));

	% Build the path/filename
	path2='d:\srp\ppt51\'; % hard coded path for data files
	pf2=[path2 c];
	disp(['Loading: ' pf2])

	% Load the .mat file  and put the monthly data (year in col 1) 
	% into A
	eval(['load ' pf2]);
	if exist('Z') & ~exist('Y1')
		A=Z;
	elseif exist('Y1') & ~exist('Z')
		A=Y1;
	elseif exist('Y1') & exist('Z')
		error([' Both Z and Y1 found in file # '  int2str(n)])
	else
		error([' Neither Z nor Y1 found in file# ' int2str(n)])
	end

	% set the start, end year of the monthl series
	yrs=YRS1(nn,:);

	% Compute the seasonalized data 
	F = seaspt(A,months(1),months(2),yrs,dtype);


	% Put seasonalized series in correct slots of Y
	i1 = F(:,1)-yr1+1;
	Y(i1,nn)=F(:,2);


end

% Compute means -- either long-term or for a uniform period;
% and sample size means based on
if O(2)==2; %  means for a uniform period
	L1=yr>=nrmpd(1) & yr<=nrmpd(2);
	nshort=sum(L1);
	Y1=Y(L1,:);
	L2=isnan(Y1);
	ns1=sum(sum(L2));
	if ns1==0; % no data has any missing data for specified period
		S(1,:)=mean(Y1);
		S(2,:)=nshort(:,ones(n,1));
	else
		for nn=1:n;
			y =Y1(:,nn);
			L3=~isnan(y);
			y=y(L3);
			S(1,nn)=mean(y);
			S(2,nn)=length(y);
		end
	end
else; % use full available data for long-term means
	for nn=1:n;
		y=Y(:,nn);
		L3=~isnan(y);
		y=y(L3);
		S(1,nn)=mean(y);
		S(2,nn)=length(y);
	end
end
		

% Normalize data, if desired
if O(1)==1'
	% no action needed, don't normalize
elseif O(1)==2;
 	s1=S(1,:);
	s2=S(2,:);
	M = s1(ones(nyrs,1),:); % means
	if dtype==1;
		Y = Y./M; % pcp as ratios, annual value to long-term mean
	elseif dtype==2; % tmp, as departures
		Y = Y - M;
	else
		error('dtype must be 1 or 2');
	end
end
		



% Organize matrices to be saved
set1=['Y yr YRS1 S months dtype O nrmpd pf1'];

[file2,path2]=uiputfile('*.mat','Store output here');
pf2=[path2 file2];
eval(['save ' pf2 ' '  set1]);


function ntot=treecnt2(pf1,t,yr1);
%ntot=treecnt2(pf1,t,<yr1>);
%
% Number of trees summed over tree-ring sites in specified
% years.
%
% Meko 3-1-97
%
%
%****************** IN  *********************************
%
% pf1(1 x :)s path,filename of file holding names of
%		sov format tree-ring storage files. Example:
% 	   c:\jack\tstfiles.txt,  holding:
%
%		c:\jack\az033
%		c:\jack\az501
%		c:\jack\dcy
%
% t (:,1)i  vector of target years 
% yr1 (1 x 1)i  optional, specifies that not only must
%		the tree exist in years specified by t, but also
%		must still be there in year yr1.  yr1 might be a 
% 		recent year (e.g., 1968), so that you want to know
%		if the tree's data could be calibrated directly 
% 		with climate data
%
%********************* OUT **************************
%
% ntot (:,1)i  number of trees in years specified by t
%
%
%***************   FUNCTIONS CALLED ******************
%
% treecnt1.m counts for individual sites
% treenum.m  analyzes string matrix of names to find unique trees
%
%*********************** NOTES ***********************
%
% Motivation: needed to count number of trees at all collection
% sites in the Galiuros in years 1600, 1650, ...


% Get the name of file with filenames
fid1=fopen(pf1,'r');

% Initial read to count number of files
k1=1;
nfiles=0;
while k1;
	c=fgetl(fid1);
	if feof(fid1);
		k1=0;
	elseif isempty(c);
		k1=0;
	elseif isspace(c(1));
		k1=0;
	else
		nfiles=nfiles+1;
	end
end


% Intialize vector to store number of trees
[nyears,nt]=size(t);	
if nt~=1;
	error('t must be cv');
end
a=NaN;
ntot=zeros(nyears,1); 

% Rewind the file of filenames
frewind(fid1);

% Loop over the files for sites
for i =1:nfiles;
	fln=fgetl(fid1);
	fln=strtok(fln); 
	% Load the mat file
	eval(['load ' fln ';']);  % will get nms and yrs
	
	% Compute number of trees for this site
	if nargin==2;  % no recent year constraint
		n=treecnt1(nms,yrs,t);
	elseif nargin==3; % recent year constraint
		n=treecnt1(nms,yrs,t,yr1);
	else
		error('Number of in args must be 2 or 3');
	end

	% Increment counter for all sites
	ntot=ntot+n;
end
fclose(fid1);

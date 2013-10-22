function [A,B]=locabn(path1,Xfn)
%
%[A,B]=locabn(path1,Xfn)
% Computes the number of locally absent (LA) rings in data file of ringwidths.
% Ringwidth file must be of a special form, which can be made using
%	dpl/yux followed by minor editing.  The required form is 
% 	row1- number of columns (variables, or cores), not counting the year
%	row2- field width of data for variables (e.g., 9), followed by headings
%		for each column. Note that the field width is above the year column
%	row3-on: the data, year in col 1, variables in other cols.  -99 as missing
% 
% CALLING FUNCTION : LOCABCML.M
% 
%*************** IN ARGS (both optional) ********************************
%
% path1 (1 x ?)s path to directory with ring-width data file
%		example:  'd:\wrk1\temp\'
% Xfn	(1 x ?)s name of file with ringwidth data, .dat assumed.
%		Example: 'bky.dat'
%
%******************* OUT AGRS ************************************
%
% A (? x 3)i  tsm telling sample size and number of LA rings, all years
%   Col 1	Year
%   Col 2	Sample size (# of cores with nonmissing data)
%   Col 3	Number of cores with ring locally absent
% B (? x 3)i	Subset of rows of A including only the rows 
%		with at least one LA value
%
%***********  USER WRITTEN FUNCTIONS NEEDED 
%
% PLTEXT1.M
% SHLMNU.M
% USINP.M
%______________________________________________________________

%  Set the path and filename of the file of ringwidth series 
if nargin==2; % path and filename are input arguments
  fln=Xfn;
elseif nargin==0; % point to the desired file
  [fln,path1]=uigetfile('*.dat','.dat file with ring widths');
else
	error('nargin must be 0 or 2');
end
pf1=[path1 fln];

% If Cancel button is pressed in the input window, flag is 
% set to -1 and control is returned to calling function
if fln==0,
  fl=-1;
  return;
end 

% open the ringwidth data file for reading
fid=fopen(pf1,'r');

% Read the total # of cores in the file from first line of data file
nstr=strtok(fgetl(fid));
ncores=str2num(nstr);
if isempty(ncores)
	fl=-1;
	disp('ncores must be a number, not a string')
 	return;
end
disp(['The total Number of columns read is : ' nstr]);


% Read the field width for data columns, then the column headings
% (variable names) from the second line of the data file 
zn=fgetl(fid);
l=length(zn);          % length of line 2
fldwd=str2num(strtok(zn));   % column width of data matrix
if isempty(fldwd),
	error('fldwd must be numeric')
end
if fldwd>15,
	error('fldwd (field width) must not exceed 15')
end


% Put the sample names in a string matrix 
bb = blanks(fldwd);
zname=bb(ones(ncores,1),:);
zn1=zn;
zn1(1:4)=[]; % lop off first 4 chars -- this would be the year field
for i = 1:ncores;
	igo=1+(i-1)*fldwd;
	isp=igo+fldwd-1;
	zname(i,:)=zn(igo:isp);
end

% Put sequence number of cores in a string matrix
nwide=length(nstr);
bb1=blanks(nwide);
num=bb1(ones(ncores,1),:);
for i=1:ncores;
	ss1=int2str(i);
	nn1=length(ss1);
	num(i,1:nn1)=ss1;
end

% Combine the name and sequence number in a string matrix
bb1=blanks(1);
zname=[zname bb1(ones(ncores,1),:) num];

% Read rest of the data file in %g format and store the data 
% in a vector
zv=fscanf(fid,'%g',inf);
tm=length(zv);
m=tm/(ncores+1);         % Number of rows in the data matrix

% Form an m x ncores matrix including the years 
zt=reshape(zv,ncores+1,m);
zt=zt';


% Build a column vector nvalid, the total number of valid (non-missing-data)
% rings in each year
zyr=zt(:,1); % cv of years
L1=zt~=-99; %  logical pointer to non-missing values in zt 
L1(:,1)=[]; %   remove col 1 of L1, which would corresp to year column
nvalid=[sum(L1')]';  % cv of sample size (number of valid values) in each yr

% Build nLA, a cv of number of LA values in each year
L2=zt==0;   % zero values (LA rings) of ringwidth
L2(:,1)=[]; % get rid of year column
nLA=[sum(L2')]'; % number of LA values in each year
A=[zyr nvalid nLA];

L3=nLA~=0;
B=[zyr(L3) nvalid(L3) nLA(L3)];

% Some years in ringwidth data file might have all -99's, meaning no valid
% data.  Want to truncate matrix A to get rid of those years.
L4=nvalid==0;
A(L4,:)=[];

% Test that no internal rows (years) of A are gone.  I.e., make sure 
% year is continuous
yr = A(:,1);
diff1=diff(yr);
L5=diff1==1;
if ~all(L5)
	error('Internal year with no valid data')
end

fclose(fid);  	% close the input file


function [A,B]=countla1(path1,file1,I,pd)
% [A,B]=countla1(path1,file1,I,pd)
%
% Count number of locally absent (LA) rings at a site in each year
%
%************************ IN ARGS **********************************
%
% path1 (1 x ?)s  path to .mat file holding ringwidth data
% file1 (1 x ?)s filename of .mat file holding ringsidth.  File assumed
%		created by rwlinp.m.  Thus, the .mat file has these variables:
%
%		X (? x 1)r the strungout vector of ringwidths
%		yrs (m x 3)i  start year, end year and row index in X of data
%			for each of m cores
%		nms(m x 8)s names of the cores
% I (mI x 1)i  index pointer to cores (rows of nms and yrs) to be used
%		in tallying the missing rings.  By default, if I is not included
%		as an argument, I==(1:m)'; this means use all cores
% pd (1 x 2)i  desired start and end year of output LA count matrix A
%		If not included, assumes that yrs covers earliest start year to
%		most recent end year
%
%
%***************************** OUT ARGS ******************************
%
% A (mA x 3)i  tally of number of valid and LA rings in each year
% 		col 1 is year
%		col 2 is number of valid (not missing data) cores
%		col 3 is number of LA values 
% B (mB x 3)i  tally as in A, but includes only those years with at
%		least one LA value
%
%
%************************* NOTES ***********************************
%
% Must run rwlinp.m previously to convert the rwl file into the .mat file


if ~exist('file1');   % Assume need to click on file
	[file1,path1]=uigetfile('*.mat','.mat file with ringwidth data');
	pf1=[path1 file1];
end

% Get the ringwidth data, core names and year info
eval(['load ' pf1]);
[mX,nX]=size(X);
[m1,n1]=size(nms);
[m2,n2]=size(yrs);

% If I not passed as input argument, assume use all m2 cores
dd=exist('I');
if dd~=1
	I = (1:m2)';
end


% Initialize matrix A
a=NaN;
if exist('pd');
	yr = (pd(1):pd(2))';
else
	yr =   ( min(yrs(:,1)):max(yrs(:,2)))';
end
nsize = length(yr);
A(:,1)=yr;
A(:,2)=zeros(nsize,1);
A(:,3)=zeros(nsize,1);


% Loop over cores
mI=length(I);
for i=1:mI;
	i1=I(i); % identifies core
	% get start and end year of this core's data
	yr1=yrs(i1,1);
	yr2=yrs(i1,2);
	igo=yrs(i1,3);
	isp = igo + (yr2-yr1);
	year=(yr1:yr2)';
	nyear=length(year);
	
	x=X(igo:isp);
	L1= year>=yr(1) & year<=yr(nsize); % will use this part of the sample
	sum1=sum(L1);
	if sum1>0
		x=x(L1);
		yrsub=year(L1);
		L2=~isnan(x); % good (not missing data?)
		L3=x==0; % LA values
		% Compute rows of A to slide this core data into
		L4=yr>=yrsub(1) & yr<=yrsub(sum1);
		% Increment valid-data counter
		A(L4,2)=A(L4,2)+L2;
		A(L4,3)=A(L4,3)+L3;
	else
	end
end


% Make output ascii files of the counts
[file2,path2]=uiputfile('*la.dat','Outfile of LA count');
pf2=[path2 file2];
fid2=fopen(pf2,'w');
fmt2='%4.0f %4.0f %4.0f\n';
fprintf(fid2,fmt2,A');








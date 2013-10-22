function [X,file1]=crn2rw(misscode)
% [X,file1]=crn2rw
% 
% Convert a Tree-Ring Lab decade-record-crn file to a ".rw"
% file that can be used by rwlook.m
%
% misscode is numeric missing value code in the crn file: like "9990"
%
[file1,path1]=uigetfile('*.crn','TR lab .crn file');
pf1=[path1 file1];
fid1=fopen(pf1,'r');


% Prelim read to get full number of lines in file
nc=0; % initialize number of lines in file to be read
k1=1;
while k1==1;
	c=fgetl(fid1);
	if ~feof(fid1);
		nc=nc+1;
	else
		k1=0;
	end
end



% Read 1st four lines to find out how many header lines.
% Consider that to be a non-header line, must have a 4-digit
% number (year) in cols 7-10; Expect 0-3 header lines.
if nc<4
	error('Fewer than 4 lines in file')
end
nh=0; % initialize number o fheader lines
frewind(fid1)
for n=1:4
	c=fgetl(fid1);
	if isempty(str2num(c(7:10)));
		nh=nh+1;
	else
	end
end
if nh>3; 
	error('More than 3 header lines in the .crn file')
end

% Rewind .crn file and skip number of header lines
frewind(fid1);
if nh>0;
	for n=1:nh;
		c=fgetl(fid1);
	end
end

% Compute adjusted number of lines to read; this many decade lines
nc=nc-nh;

% Read first data line to get the beginning year of first decade
%  e.g., 1410
c=fgetl(fid1);
decgo=str2num(c(7:9))*10;

% Compute a years vector starting with decgo and running through
% the last decade
yron=decgo;
yroff=decgo  +  (10*nc) -1;
yr1=(yron:yroff)';

% Allocate space for the year and data
a=NaN;
m1=length(yr1);
X=a(ones(m1,1),ones(2,1));
x=a(:,ones(10,1));


% Once again rewind .crn file and position after headers, if any
frewind(fid1);
if nh>0;
	for n=1:nh;
		c=fgetl(fid1);
	end
end


% Read the decade lines of data, one by one, each time converting
% the string data to numbers
%
% First set up the start and end columns for reading the 10 values
i1=[11:7:74];
i2=i1+4;
for n = 1:nc
	kgo=n*10-9;
	ksp=kgo+9;
	c=fgetl(fid1);
	% Loop over the 10 values in the decade
	for j = 1:10;
		ii1=i1(j);
		ii2=i2(j);
		x(j)=str2num(c(ii1:ii2));
	end
	X(kgo:ksp,2)=x';
end

% Put year in col 1
X(:,1)=yr1;

% Find out which leading and trailing years are missing values,
% and truncate those rows off X
L1=X(:,2)==misscode;
sum1=sum(L1);
if sum1>0;
	X(L1,:)=[];
end


% Write the 2-col ascii data file 
[file2,path2]=uiputfile('*.rw','The desired fake rw filename');
fmt3=' %g \n';
pf2=[path2 file2];
fid2=fopen(pf2,'w');

fprintf(fid2,'%s\n','fake');
fprintf(fid2,'%s\n','99/99/99');
fprintf(fid2,' %g \n',X(1,1));
xvect=round(X(:,2)/10);

fprintf(fid2,fmt3,xvect);
fprintf(fid2,' %g \n',999);



fclose all

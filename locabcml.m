function [A,B]=locabcml
%
% Cumulative sample size (number of cores) and number of LA rings
% at several sites
%
%****  IN ARGS -- NONE ********************************************** 
%
%****  IN FILES 
%
% flnms?.txt: ascii file of names of .dat files holding ringwidths for
% 	various sites; first line is path for data files, e.g., 'd:\wrk1\temp\'
%
%**************** OUT ARGS *****************************************
%
% A (? x 3)	Cumulative Locally absent indicator.
%   Col 1	Year
%   Col 2	Sample size (# of cores with valid data)
%   Col 3	Number of cores with ring locally absent
% B (? x 3)	Subset of rows of A including only the rows 
%		with at least one locally absent value
%
% USER WRITTEN FUNCTION NEEDED 
%-----------------------------
% INVTXT.M
% JDISP.M
% PLTEXT1.M
% SHLMNU.M
% USINP.M
%______________________________________________________________

[file1,path1]=uigetfile('fln*.txt','Infile with file names of .dat files')
pf1=[path1 file1]; 

% open the filename file for reading
fid=fopen(pf1,'r');

% Get the path for the data files
cpath=strtok(fgetl(fid));

b1=blanks(11);
namat=b1(ones(100,1),:); % initialize name matrix of filenames to 100 rows
nfiles=0;
while 1
	% get rid of leading and trailing blanks
	c=fgetl(fid);
	if feof(fid);
		break
	end
	nfiles=nfiles+1;
	c=strtok(c);
	c=fliplr(c);
	c=strtok(c);
	c=fliplr(c);
	clen=length(c);
	namat(nfiles,1:clen)=c;
end
namat=namat(1:nfiles,:);
disp(namat)
fclose(fid);

[na,nb]=size(namat); % na is number of files, nb is col size of names matx

rnl=zeros(na,1);
yrbeg=zeros(na,1);
yrend=zeros(na,1);
for i=1:na,
  rnl(i)=length(namat(i,:));
  yrbeg(i)=str2num(namat(i,rnl(i)-13:rnl(i)-9));
  yrend(i)=str2num(namat(i,rnl(i)-4:rnl(i)));
end

ybg=min(yrbeg);
ynd=max(yrend);
A=zeros(ynd-ybg+1,3);  		% Preallocate memory for A and B
B=zeros(ynd-ybg+1,3);
A(:,1)=(ybg:ynd)';
B(:,1)=(ybg:ynd)';

for i=1:na,
  ofnam=namat(i,1:16);
  L1=ofnam'==' ';
  ofnam(L1)=[];
  [AT,BT]=locabn(ofnam);	% Find A and B matrices for each file
  lnat=length(AT(:,1));

  yr1=AT(1,1);
  yr2=AT(lnat,1);
  lga=A(:,1)>=yr1 & A(:,1)<=yr2;
  A(lga,2:3)=A(lga,2:3)+AT(:,2:3);	% Cumulative sum of 2nd & 3rd columns

end

la=find(A(:,2)==0);
A(la,:)=[];			% Take out all the zeros from A and B
B=A;
lb=find(B(:,3)==0);
B(lb,:)=[];






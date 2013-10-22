function chkyrs1
%
% Given (1) a group of monthly pcp  *.mat files,  
% (2) a txt file of names of the .mat files, and (3) a .mat file
% with supposed start and end years for the .mat files --- does this:
%
% Checks that the indicated start and end years match the
% true start and end years in the files. 
%
% Meko 8-13-96
%
% Written to aid in verifying the pcp.dbf database info on 
% time coverage of estimated pcp data for SRP study
%
%********************* ASSUMES
%
% The monthly data .mat files are in the current directory
% The file of file names and of supposed first,last years is too
% The order of stations in the filename and years files is the same

clear

% Load the .mat file of start, end years
[file2,path2]=uigetfile('year*.dat','File with start, end years');
pf2=[path2 file2];
eval(['load '  file2]);
file2=strtok(file2,'.');
eval(['YRS = ' file2]);
[m1,n1]=size(YRS); % n1 is number of stations

 %Open the file of file names
[file1,path1]=uigetfile('*.txt','File with .mat file names');
pf1=[path1 file1];
fid1=fopen(pf1,'r');


n=0; % initialize counter
k1=1;
while k1==1;
	c=fgetl(fid1);
	if ~feof(fid1)
		n=n+1;
		yr1=YRS(n,1);
		yr2=YRS(n,2);
		c=strtok(c);
		c=c(~isspace(c));
		% load the file
		eval(['load ' c])
		if exist('Y1') & ~exist('Z')
			Z=Y1;
			disp('a file with Y1')
			eval(['save ' c ' Z'])
		end
		[m2,n2]=size(Z);
		yra=Z(1,1);
		yrb=Z(m2,1);
		if yr1==yra & yr2==yrb;
			disp(['Series ' int2str(n) ' years OK']);
		else
			disp(['Series ' int2str(n) ' years NOT OK']);
			c
			x=[yr1 yr2 yra yrb];
			disp(x)
			pause
		end
	else
		k1=0;
	end
	clear Z
	clear Y1
end
		


fclose (fid1)















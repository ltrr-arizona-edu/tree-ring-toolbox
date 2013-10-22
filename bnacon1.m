function bnacon1
%
% Cut down a mapview .bna file to include only points in a 
% specified long-lat box
%
% Meko 12-6-96
%
%***************** IN FILES *********************************
%
% 1. Ascii file with two lines.  Bounding lats on line one, longs
% on line two. Sample:
%
% 38.8 40.3
% 122.2 124.6
%
% 2. ..bna file as produced by export from mapviewer
%
%
%**************** OUT FILE
%
% Ascii .bna file, trimmed to size


% Get the lat, long specs
[file1,path1]=uigetfile('*.dat','Input lat/lon limits file');
pf1=[path1 file1];
eval(['load ' pf1]);
eval(['x = ' strtok(file1,'.') ';']);
lat1=x(1,1);
lat2=x(1,2);
lon1=-x(2,1);
lon2=-x(2,2);

% Open the .bna file
[file2,path2]=uigetfile('*.bna','Input bna file');
pf2=[path2 file2];
fid2=fopen(pf2,'r');

% Read thru the input bna file, building a logical vector same size
% as number of lines in it, with 1 meaning this is a bna header line, 0
% a data point
%
% Allocate
L1=ones(50000,1);  % hard code max number of lines
nlines=0;
nc = 0;
k1=1
while k1
	c=fgetl(fid2);
	if feof(fid2)
		% no increment to nlines
		k1=0;
	else
		nc=nc+1;
		if c(1)=='"';
			L1(nc)=1;
		else
			L1(nc)=0;
		end
	end
end; % of while k1	
L1	=L1(1:nc); % truncate unused elements of L1

% Find row number of each bna header in L1
ihead = find(L1);
nhead = length(ihead);

% Get start and end row number for x,y points in each header section
igo=ihead+1;
ih2=ihead;
ih2(1)=[];
isp=[ih2-1; nc];
ntot=isp-igo+1; % number of x,y points in each section


% Open file for output bna
[file3,path3]=uiputfile('*.bna','Output bna file');
pf3=[path3 file3];
fid3=fopen(pf3,'w');


%********* LOOP OVER HEADER SECTIONS, reading the data points and writing
% to a new file only those in the desired long-lat range. Also keep track
% of how many culled points in each section. If fewer than two points, do
% not write for that section

frewind(fid2); % rewind the bna file

for n = 1:nhead;
	disp(['On ' int2str(n) ' of ' int2str(nhead)]);
	np=0; % counter for number of valid points
	n1=ntot(n); % number of points before culling
	c1=fgetl(fid2); % should get the header line
	S1=c1;
	for m = 1:n1;
		c2=fgetl(fid2);  % sting with long, comma, lat
		x=str2num(strtok(c2,',')); % decimal neg long
		fcomma=findstr(c2,',');
		c3=c2;
		c3(1:fcomma)=[];
		y=str2num(strtok(c3)); % decimal deg lat
		if x>=lon2 & x<=lon1 & y>= lat1 & y<=lat2;
			S1=str2mat(S1,c2);
			np=np+1;
		end
	end
	% If two or more points in S1, write it to file
	if abs(np)>1;
		% Change header line specification of number of points, if needed,
		% and to negative so that surfer will not expect a closed area
		if abs(np)<abs(n1); % some points have been deleted
			if np>0; % used to be an area, now change to a curve
				np1=-1*np;
			end
			hd = S1(1,:);
			msize=length(hd);  % is also col size of S1
			f2 = findstr(hd,',');
			sbeg=hd(1:max(f2));  % part of header from start thru last comma
			strnp=num2str(np1); % string version of number of valid points
			len1=length(strnp);
			len2=msize-length(sbeg)-len1; % will need to left-pad this many blanks
			S1(1,:)=[sbeg blanks(len2) strnp];
		end
		fprintf(fid3,'%s\n',S1(1,:));
		S1(1,:)=[];
		for j = 1:np;
			s1=S1(j,:);
			fprintf(fid3,'%s\n',s1);
		end
	else; % fewer than 2 points, no action needed
	end
end		



fclose all



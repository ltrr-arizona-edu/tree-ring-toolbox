function [nstns,S,N]=tmpcon2
% tmpcon2: convert a NCDC cooperative stations file of monthly ppt series
% to individual station .mat files;  also build string files of
% station names and filenames, and a numeric file of first and last
% year; also build ascii tab-delimited file for station database
% CALL: function [nstns,S,N]=tmpcon2;
%
% D Meko 6-16-97 -- modeled on the corresp ppt function, pptcon2.m
%
% Assumed form of input file -- ftp'd from NCDC
%
% * 1-line header with "NUMBER" in cols 2-7
% * data lines, one after another
%	ch 2-7 station id (e.g., 020108 for an AZ station) 
% ch 10-32 station name
%	ch 33-34 state (AZ or NM for SRP study)
%	ch 40-43 year
%	ch 44-49, 50-55,... 110-115 etc -- jan ... dec. pcp 
%		in hundredths of inches
%	ch 116-123  annual total 
%
% * 9999 right justified in field means missing data
% * possible that entire internal years might be missing
% * '$' in col 1 of line after last line of data to be read in
%
% Assumed output form
%
% * Data
% 13 cols :  year followed by 12 values in inches, with decimal point
% Missing data as NaN
% Any internal years of missing data filled with NaN rows
%
% * Data filenames
% Prefix from station id, with AZ replacing "02" and NM replacing
%	"29".  Root ended with "p" for pcp.  
%  Example AZ0040p.mat for 020040, Agua Caliente, AZ
%
% * filename-list fil -- ? x 12, right-filled with blanks
% * station name   -- (? x 23)s
% * start,end  file -- (? x 2)i. first, last year
%
% DATABASE FILE -- see comments in hcndbf1.m for contents of
%   12 columns


%**************** get input file

% Work with Arizona or New Mexico?
k2=input('Arizona (1) or New Mexico (2) data:  ');
if k2==1,
	ste='AZ';
elseif k2==2
	ste='NM';
else
	error('k2 must be 1 or 2');
end

% Build name of input data file
if k2==1; % AZ
	file1='arizona.txt';
else
	file1='new_mexi.txt';
end

% Click on and open the input data file
[file1]=uigetfile(file1,'Input Coop monthly NCDC file');
fid1=fopen(file1,'r')
frewind(fid1)

% Make temporary file for intermed output
file2='temp001.dat';
fid2=fopen(file2,'w')


a=NaN;


maxn=1000; % for now, allow for max of 1000 stations

igo=a(ones(maxn,1),:); % row index of starting year for station

%Initialize matrix of station names
s3=blanks(23);
N = s3(ones(maxn,1),:);

% Initialize string matrix of file names
s1=blanks(12);
S=s1(ones(maxn,1),:); 

% Initialize matrix for start end years
YRS=a(ones(maxn,1),ones(2,1));


% Read the 1-line header
c=fgetl(fid1);
if all(c(2:7)=='NUMBER');
	disp('One line header read')
else
	error('Header line should have NUMBER in cols 2-7')
end

% Strategy. A new series is being read if the id in cols 2-7 for the
% line being read is different from the previous line
% -initialize idprev
% -read a line and compare id to idprev
% -if same, increment line counter for this station
% -if different, you have finished a series
%	-increment series counter
%	-set ending row as last line # read
%	-set ending year as year of prev line read
%	-move back a line in the file
%	-decrement the line counter
%	-set idprev==id
%	-set start year and row


knew=1; % ==1 would mean have just read a line 1
kln2=0; % ==1 would mean have just read line 2
k1=1; % control for while loop

lineno=0; % initialize line counter
stnprev='12345678901234567890123';
idprev='000000';


nstns=0; % initialize number of stations
while k1==1
	c = fgetl(fid1);
	%disp(c)
	if c(1)~='$'; % flag for end of file
	%if c(1)~='$' & nstns<4; % flag for end of file
		id=c(2:7);
		stn=c(10:32); % station name
		cdata=c(40:115);  % all-numeric part, from year through dec value
		if lineno==0 |  ~all(id==idprev) ;  % starting a new station
			knew=1; % starting a new station
			j=0; 
			nstns=nstns+1;
			str13=sprintf('  Station no. %5.0f',nstns);
			cfile=[ste strtok(c(4:7),' ') ]; % will be the .mat file
			disp(['Part 1, working on : ',cfile,str13]);
			nss = length(cfile);
			S(nstns,1:nss)=cfile;
			N(nstns,:)=stn;
			lineno=lineno+1;
			igo(nstns)=lineno;
			fprintf(fid2,'%s\n',cdata);
		else; % another line for this station
			lineno=lineno+1;
			fprintf(fid2,'%s\n',cdata);
		end
		idprev=id;
	else; % hit the "$"
		k1=0;
	end
end;  % of while k1 block

fclose(fid1);
fclose(fid2);

ntot=lineno;

istop=igo; % recall that igo is row index of start year
igo=igo(1:nstns);
istop=istop(1:nstns);


istop(1)=[];
istop=istop-1;
istop=[istop; ntot]; % row index of ending year 

% Lop off unneeded rows
S = S(1:nstns,:);
N = N(1:nstns,:);

clc
disp(['Total of ' int2str(nstns) ' stations']);
disp(['Total of ' int2str(ntot) ' lines read']);


%*******************************************************************
%
% First part finished.  temp001.dat is an intermediate file without
% the header lines, and with only the numeric data of year and
% 12 monthly values
%
% temp001 still has 9999 as missing, and might have entire years 
% missing

%**********************************************
disp('Starting second part')
load temp001.dat;
X=temp001;

% Loop over stations
for n = 1:nstns;
	stnname=N(n,:);
	disp(['working on ' stnname]);
	Y=X(igo(n):istop(n),:);  % data for a station
	yr =Y(:,1); % year vector
	A=Y(:,2:13); % pcp in hundredths of inches

	% replace 9999 with NaN
	L1=A==9999;
	sum1=sum(sum(L1));
	if sum1>0;
		A(L1)=a(ones(sum1,1),:);
	end

	% Scale data from hundredths of inches to inches
	A=A* 0.01;
	

	yrgo =min(yr); % lowest year for input
	yrsp=max(yr); % highest year for input
	yrr = (yrgo:yrsp)';  % year vector for output
	nyears = length(yrr);  % Number of years in output matrix
	YRS(n,:)=[min(yrr)  max(yrr)];
	Z=a(ones(nyears,1),ones(12,1));

	% Fill matrix Z with the scaled and NaN-filled data
	rowi = yr - yrgo +1;  % row index in Z for data in A
	Z(rowi,:)= A;

	% put back the years column
	Z = [yrr Z];

	% Build output filename
	fnout = [strtok(S(n,:)) 'p'];
	
	% Save the .mat file
	eval(['save ',fnout,' Z']);
end			
						

%*********************************************************
% Build site info database ascii file    ??cdbf.txt, where
% ?? is az or nm
disp('Starting third part: building site info file')

% coop2.txt is assumed to hold the needed site info 
fid3=fopen('coop2.txt','r');

% azdbf.txt or nmdbf.txt will hold the output ascii dbf info
if all(ste=='AZ');
	file4='azdbf.txt';
elseif all(ste=='NM');
	file4='nmdbf.txt';
else
	error('ste must be AZ or NM')
end
fid4=fopen(file4,'w');

str3=ste;
str11='c'; % a "cooperative" station

fmt4='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';

% Loop over stations
for n = 1:nstns;
	str12=sprintf('%3.0f',n);
	%disp(['n = ',int2str(n)]);
	%disp(['str12 = ' str12]);
	str1=[deblank(S(n,:)) 'P'];;  %  file name, without .mat
	% Get 2-number state code
	if all(ste=='AZ');
		dd='02';
	elseif all(ste=='NM')'
		dd='29';
	else
		error('ste must be 02 or 29')
	end


	% Find the right station info in coop2.txt
	this=[dd S(n,3:6)];   
	k7=1;
	while k7==1;
		c=fgetl(fid3);
		if(feof(fid3));
			fclose(fid3);
			fclose(fid4);
			disp(['While on series ' this]);
			error('Reached end of file on fid3'); 
		end
		%disp(['Station number ' int2str(n)]);
		%disp(c)

		cstate=c(22:23);
		cfour=c(24:27);
		cdiv=c(29:30);
		that = c(22:27);
		%disp(['this = ' this])
		%disp(['that = ' that])
		if all(this==that); % found first id entry matching key
			% Lat and long to decimal plotting degrees
			disp(['Found match for series ' this]);			

			str2=c(1:21);

			clatdeg=c(50:51);
			clatmin=c(52:53);
			c1=str2num(clatdeg);
			c2=str2num(clatmin);
			c3=c1+c2/60;
			c4=(round(c3*100))/100;
			clat=sprintf('%5.2f',c4);

			clondeg=c(56:58);
			clonmin=c(59:60);
			c1=str2num(clondeg);
			c2=str2num(clonmin);
			c3=c1+c2/60;
			c4=(round(c3*100))/100;
			c4=-1.0*c4;
			clon=sprintf('%6.2f',c4);


			str4=cfour;
			str5=cdiv;
			str6=clon;
			str7=clat;

			celev=c(62:67); % elev in tenths of meters
			c1=str2num(celev);
			c2=(c1/10); % elev in meters
			str8=sprintf('%7.1f',c2);

			str9=sprintf('%4.0f',YRS(n,1));
			str10=sprintf('%4.0f',YRS(n,2));
			fprintf(fid4,fmt4,str1,str2,str3,str4,str5,...
				str6,str7,str8,str9,str10,str11,str12);
			k7=0;  % to get out of while loop
		else; % not a match in id's -- do nothing, stay in while loop
		end
	end; % of while loop

end; % of for n=
	
fclose(fid3);
fclose(fid4);


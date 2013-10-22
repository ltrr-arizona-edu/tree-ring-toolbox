function [nstns,S,P]=pptcon3
%
% read a lofgren-style file of multiple monthly ppt series;
% convert to individual station.mat files of monthly data;
% convert units from hundredths of inches to inches
% build a site-info database specific for SRP study
% build other misc numeric files
%   - long/lat, elev (m) plotting file
%	 - station filenames file
%
% Modeled on pptcon1.m
%
%  D Meko  7-15-96
%
%
% Assumed form of input file
%
% * sequential station series with two header lines
%   header-1:  '*' in col 1, followup by desired file
%			name of up to 5 chars, "." and suffix not included
%	 header-2:
%			ch 1-2    state (02 for arizona, etc)
%			ch 3-4    clim division number
%			ch 5-8    4-character station id
%			cy 13-39  station name
%			ch 44-46  'PPT'
%			ch 66-67  lat deg
%			ch 68-69  lat min
%			ch 71-73  lon deg
%			ch 74-75  lon min
%			ch 77-80  elev (ft)
% * data lines, one after another
%	ch 1-4	year, 
%  ch 5-64 monthly values in A5, hundreths of inches
%	ch 65+ trailing garbage no to be used
%
% * blank 5 slots means missing monthly data
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
% Prefix from leading header line, plus "b", suffix .mat
% where the "b" indicates the "Bradley" data set
%  Example *tucs -->  tucsb.mat
%			*ajo --> ajob.mat
%
% * filename-list fil -- ? x 8 right-filled with blanks
% * database info file
% Want a new file braddbf.txt with 10 cols that will be fields in 
%   site info database;  will match fields created for
%	 cooperative stations by pptcon2.m, and for USHCN stations by
%	 hcndbf1.m
%
% Want this order of cols
%
% 1 file name, with "b", without suffix  --e.g., ajob
% 2 station name
% 3 state -- AZ or NM
% 4 id -- 4-character, all numbers
% 5 div -- *** this will be set at dummy "99" because hcn file
%   info does not include clim division.  Will need to put in 
%   database later by hand
% 6 lon -- decimal deg, neg w
% 7 lat -- decim deg
% 8 el_m  elevation in meters
% 9 yrgo1 - first year with data
% 10 yrsp1 - last year with data
% 11 src -- source of data (=='b') "Bradley"
% 12 seq -- sequential number of station in Bradley save set ppt.bra


%**************** get input file

file1=uigetfile('ppt.bra','Input Lofgren-style PPT file Bradly data');
fid1=fopen(file1,'r')
frewind(fid1)

file2='temp001.dat';
fid2=fopen(file2,'w')

a=NaN;

% Initialize string matrix of file names
s1=blanks(8);
s2=blanks(27);
s3='  '; % state
s4='    '; % id
s5='  '; % division
S=s1(ones(200,1),:); % assume max of 200 stations treated
S2=s2(ones(200,1),:); % station names
S3=s3(ones(200,1),:); %state
S4=s4(ones(200,1),:); % id
S5=s5(ones(200,1),:); % division

YRS=a(ones(200,1),ones(2,1)); % for start and end year

P=zeros(200,3); % to hold numeric lon,lat,elev (m)
LL=zeros(200,2); % start and end line numbers of data for each
		% series in intermediate file

knew=1; % ==1 would mean have just read a line 1
kln2=0; % ==1 would mean have just read line 2
k1=1; % control for while loop

lineno=0;

nstns=0; % initialize number of stations
while k1==1
	c = fgetl(fid1);
	if c(1)=='*';  % starting a new station
		knew=1; % starting a new station
		j=0; 
		nstns=nstns+1;
		c(1)=[]; % lop off the '*'
		ss = strtok(c,' ');
		disp(['Part 1, working on : ',ss])
		nss = length(ss);
		S(nstns,1:nss)=ss;
		% Initialize the data storage matrix
			X=a(ones(200,1),ones(13,1));
	elseif c(1)=='$'; % marker after last data line
		k1=0;
	else
		if length(c)<44;
			disp(['Station ' ss]);
			error('Not a * or $ line, but length <44')
		end
		if isletter(c(44)); % must be a "line 2"
			S2(nstns,:)=c(13:39); % station name


			cstate=c(1:2); % state code
			if all(cstate=='02');
				S3(nstns,:)='AZ';
			elseif all(cstate=='29');
				S3(nstns,:)='NM';
			else
				error('should only have AZ and NM')
			end
			

			S4(nstns,:)=c(5:8);  % station id
			S5(nstns,:)=c(3:4); % clim division

			latdeg=str2num(c(66:67));
			latmin=str2num(c(68:69));
			londeg=str2num(c(71:73));
			lonmin=str2num(c(74:75));

			% Decimal plotting long and lat
			P(nstns,1)=  -1.0* (londeg+lonmin/60);
			P(nstns,2)=  latdeg +latmin/60;

			% Input elev is in ft, want output in meters
			xelev=str2num(c(77:80));
			P(nstns,3) = str2num(c(77:80))/3.2808;
		else; % must be a data line
			lineno=lineno+1;
			if knew==1
				LL(nstns,1)=lineno;
				knew=0;
			end
			cdat=c(1:64);
			fprintf(fid2,'%64s\n',cdat);
		end
	end
end;  % of while k1 block
fclose(fid1);
fclose(fid2);
ntot=lineno;
nlast=lineno+1;
m1=[LL(1:nstns,1); nlast];
nlines= diff(m1); % number of lines for each station in file2


S = S(1:nstns,:);
P = P(1:nstns,:);
LL=LL(1:nstns,:);
LL(:,2)=LL(:,1) + (nlines-1);

%*******************************************************************
%
% First part finished.  temp001.dat is an intermediate file without
% the header lines, and with only the first 64 columns of the orig
% data file. But temp001.dat is still not guaranteed numeric loadable
% because no space between monthly values in some series.  Need to
% convert.  Also need to handle blank entries and make the NaN


disp('Starting second part')

X=a(ones(ntot,1),ones(13,1));
fid1 = fopen(file2,'r');

x=zeros(1,12);
lineno=0;

% Loop over stations
for n = 1:nstns;
	disp(['Station # ' int2str(n)]);
	% Loop over years
	for k = 1:nlines(n);
		c = fgetl(fid1);
		lineno=lineno+1;
		yr = str2num(c(1:4));
		for j=1:12;
			jgo = (j-1)*5+5;
			jsp = jgo+4;
			xtemp=str2num(c(jgo:jsp));
			if isempty(xtemp),
				x(j)=NaN;
			else
				x(j)=xtemp;
			end
		end
		X(lineno,:)=[yr x];
	end
end			
						

%****************************************************************
%
% Second part done.  Now have all-numeric matrix X, with
% year in col 1, 12 monthly values in other columns.  But 
% what if an entire year's or several years' data are missing.
% Must make NaN rows then

disp('Starting part 3')

% Loop over stations
for n = 1:nstns
	igo=LL(n,1);
	isp=LL(n,2);
	Y = X(igo:isp,:);
	Y= [Y(:,1)  Y(:,2:13)/100]; % convert to units of inches
	yr = Y(:,1);
	yr2=(min(yr):max(yr))';
	YRS(n,:)=[min(yr) max(yr)];
	nyears = length(yr2);

	if nyears ~= (isp-igo+1); % an internal year missing?
		disp(['Station ',S(n,:)]);
		disp('Computed that an internal year must be missing')
		pause(2);
		Z=a(ones(nyears,1),ones(13,1));
		i1 = yr-yr2(1)+1;
		Z(i1,:)=Y;
	else
		Z=Y;
	end

	% Build output filename
	fnout = [strtok(S(n,:)) 'b'];
	
	% Save the .mat file
	eval(['save ',fnout,' Z']);
end
		


%*********************************************************
% Build the site-info database string file

disp('Starting part 4: build database file');



file5=uiputfile('bradsrp.dbf','Output .dbf file');
fid5=fopen(file5,'w');
fmt5='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';

str11='b'; % source is Bradley data set

for n=1:nstns;
	str12=sprintf('%3.0f',n); % sequential number
	% Build output filename
	str1 = [strtok(S(n,:)) 'B']; % name of .mat data file
	str2=S2(n,:); % station name
	str3=S3(n,:); % state abbrev
	str4=S4(n,:); % station id
	str5=S5(n,:); % division
	
	str6=sprintf('%7.2f',P(n,1)); % long
	str7=sprintf('%7.2f',P(n,2)); % lat
	str8=sprintf('%7.1f',P(n,3)); % elev in m

	str9=sprintf('%4.0f',YRS(n,1)); % start year
	str10=sprintf('%4.0f',YRS(n,2)); % end year

	fprintf(fid5,fmt5,str1,str2,str3,str4,str5,str6,...
 		str7,str8,str9,str10,str11,str12);
end

fclose (fid5);

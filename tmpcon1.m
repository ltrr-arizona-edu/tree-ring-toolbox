function [nstns,S,P]=tmpcon1
%
% convert a lofgren-style file of multiple monthly hcn tmp series
% to individual station .mat files;  also build a string file of
% station-file names, string file of ids, and numeric file of
% lat,long,elevation (ft)
%
% Assumed form of input file
%
% * sequential station series with two header lines
%   header-1:  '*' in col 1, followup by desired file
%			name of up to 5 chars, "." and suffix not included
%	 header-2:
%			ch 1-8:  station id
%			ch 49-52  'Temp'
%			ch 66-67  lat deg
%			ch 68-69  lat min
%			ch 70-72  lon deg
%			ch 73-74  lon min
%			ch 77-80  elev (ft)
% * data lines, one after another
%	ch 1-4	year, 
% ch 5-64 monthly values in A5, hundreths of deg F
%	ch 65+ trailing garbage no to be used
%
% * user prompted for missing value code <-9999>
% * blank 5 slots assumed missingvalues
% * possible that entire internal years might be missing
% * '$' in col 1 of line after last line of data to be read in
%
% Assumed output form
%
% * Data
% 13 cols :  year followed by 12 values in deg F, with decimal point
% Missing data as NaN
% Any internal years of missing data filled with NaN rows
%
% * Data filenames
% Prefix from leading header line, plus "t", suffix .mat
%  Example *tucs -->  tucst.mat
%			*ajo --> ajot.mat
%
% * filename-list fil -- ? x 12, right-filled with blanks
% * station ID file   -- ? x 8
% * location file -- (? x 3)i;  decimal degrees, plotting format
%
%



%**************** get input file

file1=uigetfile('ushcn.mad','Input Lofgren-style HCN tmp file');
fid1=fopen(file1,'r');
frewind(fid1)

% Get missing value code
prompt={'Enter code'};
def={-9999};
title='Missing Value Code in Input';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
msscode=str2num(answer{1}); 	

% Scale from hundredths of deg F to deg F
prompt={'Multiply by this: '};
def={0.01};
title='Scaling factor to multiply monthly values by';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
fact1=str2num(answer{1});



file2='temp001.dat';
fid2=fopen(file2,'w');


a=NaN;

% Initialize string matrix of file names
s1=blanks(12);
S=s1(ones(200,1),:); 

P=zeros(200,3);
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
		if length(c)<41;
			disp(['Station ' ss]);
			error('Not a * or $ line, but length <40')
		end
      if isletter(c(41)); % header line 2 or 3
         if c(41)=='M' ; % must be a "line 2"
            latdeg=str2num(c(66:67));
            latmin=str2num(c(68:69));
            londeg=str2num(c(70:72));
            lonmin=str2num(c(73:74));
            P(nstns,3) = str2num(c(77:80));
            P(nstns,1)=  -1.0* (londeg+lonmin/60);
            P(nstns,2)=  latdeg +latmin/60;
         elseif c(41)=='w'; % must be line 3 header
            % no action needed
         else
            error('Invalid col 40 entry on header line 2 or 3');
         end
         
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
% convert.  Also need to handle blank entries and missing value codes
% and make them NaN


disp('Starting second part')

X=a(ones(ntot,1),ones(13,1));
fid1 = fopen(file2,'r');

x=zeros(1,12);
lineno=0;

% Loop over stations
for n = 1:nstns;
   disp(['   Working on ' S(n,:)]); 
	% Loop over years
	for k = 1:nlines(n);
		c = fgetl(fid1);
		lineno=lineno+1;
		yr = str2num(c(1:4));
		for j=1:12;
			jgo = (j-1)*5+5;
			jsp = jgo+4;
			xtemp=str2num(c(jgo:jsp));
			if isempty(xtemp) | xtemp==msscode;
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
	yr = Y(:,1);
	yr2=(min(yr):max(yr))';
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
   
   
  Z(:,2:13)=Z(:,2:13)*fact1;
  
	% Build output filename
	fnout = [strtok(S(n,:)) 't'];
	
	% Save the .mat file
	eval(['save ',fnout,' Z']);
end
		

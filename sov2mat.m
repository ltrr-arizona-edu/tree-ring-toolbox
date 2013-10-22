function sov2mat
% multiple time series in strung-out-vector for to time series
% matrix
%
% D Meko 4-22-96
%
%
%************** INPUT -- SCREEN PROMPTED **********************
%
% The user will be prompted for file names for the following input:
%
% file1 -- .dat: V (mV x 1)r:  the strung-out vector
% file2 -- .dat: P (mP x 5)i: pointer matrix to years
%		and rows for each time series in V. Cols of P are
%		1-sequence no. 2,3-start,end year, 4,5-start,end row index
% 
% file3 -- .txt: T(mT x nT)s, the name-tag
%		file from which series id or names can be drawn. 
%		First 6 characters of these names are stored in string matrix
%		along with the time series matrix in the ouput .mat file
%		These names can also be used as column headings ifthe user
%		wants a BMDP-type output matrix
% file4 -- .dat: I (mI x 1)i row pointer to T indicating the subset
%		of series to be included in the output matrix
% file5 -- name of .mat file to hold the output time series
%		and its associated names of series 
% file6 -- (optional) name of output .dat file to hold the BMDP-style
%		output matrix -- not yet built in
%
%  OTHER SCREEN-PROMPTED INPUT
%
% yrgo,yrsp -- desired start and end year of the output matrix
% xmiss -- filler code for missing values in the output matrix
%		Not yet built in; program now uses NaN
%
%************** OUTPUT **********************
%
% 1 -  .MAT  file holding matrices X and N.  X is the time-series
%		matrix for years yrgo to yrsp with data for the subset
%		of variables identified by I. First row of X is the year.
%		N is the string matrix of series names for X.
% 2 (optional) -  .DAT holding same information as the .MAT file,
%		but with the names as column headings
%
%****************** NOTES ***************************************
%
% Leading 8 chars of the lines in the name-tag file are used
% as the series identifiers.  If fewer than 8 characters 
% (e.g.,   x01.crn), right-padded with blanks




%
file1=uigetfile('*.dat','strgout_ vector input data file');
file2=uigetfile('*.dat','rowyr_ pointer matrix for input sov');
file3=uigetfile('*.dat','crnfn_ chronology file names input file');
file4=uigetfile('*.dat','use_  series-to-use pointer file');


% Load the strung-out-vector file
eval(['load ',file1]);

% Extract the matrix name for the strung-out vector; store the 
% strung-out vector as V
line=file1;
i1=find(line=='.');
vtemp=line(1:(i1-1));
eval(['V = ',vtemp,';'])
eval(['clear ',vtemp]);


% Load file with the years/rows pointer
eval(['load ',file2]);

% Extract the matrix name for the years pointer; store in P 
line=file2;
i1=find(line=='.');
vtemp=line(1:(i1-1));
eval(['P = ',vtemp,';'])
eval(['clear ',vtemp]);

% Check P
[mP,nP]=size(P); % mP is the number of time series in the sov
if nP~=5,
	error('P should have 5 columns')
end


% Get the series name or id information off a .txt file.  
% This file holds a string in each line.  The string identifies
% a time series in the sov.  The strings might differ in length
% from line to line, so cannot read in a matrix.
fid3=fopen(file3,'r+');
n1=0; % counter for number of lines in file3 
while 1
	line = fgetl(fid3);
	if ~isstr(line), break, end
	n1=n1+1;
end

if n1~=mP,
	error('N of lines in file3 not equal to row-size of P')
end

fseek(fid3,0,-1); % rewind file3

% Build name-tag string matrix S.  Assume col-size 8
S='        ';
for n = 1:n1;
	line=fgetl(fid3);
	line=deblank(line);
	ns = length(line);
	if ns>8,
		line=line(1:8);
	end
	S=str2mat(S,line);
end
S(1,:)=[];
fclose(fid3)

% Load the file with the pointer to rows of P and S indicating
% which series to put in the output matrix
eval(['load ',file4]);

% Extract the matrix of the selection pointer; store in I 
line=file4;
i1=find(line=='.');
vtemp=line(1:(i1-1));
eval(['K= ',vtemp,';'])
eval(['clear ',vtemp]);
[mK,nK]=size(K);

if nK~=1,
	error('File 4 should consist of single column')
end
if mK>mP 
	error('Row-size of K too large')
end
if min(K)<1 | max(K)>mP,
	error('A value of K in file4 is too big')
end


% Cull the desired rows of S,P
N=S(K,:); % names of desired series
P=P(K,:);

% What time period do you want the output to cover?
igo=input('Start year of output matrix: ');
isp=input('End year of output matrix: ');
if isp<=igo,
	error('Stop year must be later than start year')
end
if igo>=max(P(:,3))| isp <= min(P(:,2)),
	error('isp,igo inconsistent with time coverage in P')
end
nyrs = isp-igo+1;


% Initialize a matrix to store the output tsm 
a = NaN;
X=a(ones(nyrs,1),ones(mK+1,1));
yrX = (igo:isp)';
X(:,1)=yrX;

% Loop over the desired series, pulling out the rows (years)
% in the "igo to isp" period, if any, and putting in the correct
% rows of X .  Must allow for possibility that some selected
% series may have no non-missing data in the nyrs key period.
for n = 1:mK;
	yr1 = [P(n,2):P(n,3)]';
	nsize=P(n,3)-P(n,2)+1; % number of years of data for this series
	if (yr1(1)>isp | yr1(nsize)<igo); % No usable data 
		% skip this series
	else
		v = V(P(n,4):P(n,5)); % cv if data in V for this series
		% Make two pointer logical vectors.  L1 is pointer to v
		% marking rows in the igo to isp target period. L2 is a
		% pointer to X marking the usable years from V for this
		% series
		L1 = yr1>=igo  & yr1<=isp; % pointer to desired rows of v
		L2 = yrX>=P(n,2) & yrX<=P(n,3); % 
		% Put the data in the slots
		X(L2,n+1)=v(L1);
	end
end
%

% Save result in .mat file
file5=uiputfile('*.mat','Save tsm and names string matrix as: ');
eval(['save ',file5,' X N'])



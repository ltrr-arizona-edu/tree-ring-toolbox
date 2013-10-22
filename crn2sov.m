function [v,YRS] = crn2sov
% crn2sov: multiple tree-ring series from .crn files to a single strung-out vector
% CALL: crn2sov
%
%*********************  IN **************************%
%
% No input arguments, but the following files required
%
%		file1 - input, control file, defined by example below:
% 		669............ number of .crn files to read
% 		1850 1995 ..... include data from this year range only in v (see OUT)  
% 		c:\wrk6\crnflsa.dat ....  ascii file with path\filename for each .crn file
% 		c:\wrk6\gospa.dat ...start and end year for each chronology
% 		c:\wrk6\strgouta.dat...outfile <optional> ascii version of v
% 		c:\wrk6\rowyra.dat...<optional> outfile of year and row reference in v
% 		c:\wrk6\dummy.dat
%
%		file2 - input, names of .crn files -- from foxpro (crnflsa.dat in example above)
%		file3 - input, start, end years of chronology coverage (gospa.dat in example above)
%		.crn files:  many files sequentially opened, operated on, and closed
%
%********************* OUT ***************************
%
% v (? x 1)r strung-out vector
%		file6 - output, strung-out vector
% YRS (nsers x 5)i year and row in sov for nsers chronologies
%		col 1: sequential number of file
%		col 2: start year of series
%		col 3: end year of series
%		col 4: row in sov of first value for this series
%		col 5: row in sov of last value for this series
%
%*********************  NOTES ************************************
%
%
% INCLUDE FULL PATHS IN FILESNAMES IN FILE2 IF ALL 
% CHRONOLOGY FILES ARE NOT IN THE DEFAULT CURRENT DIRECTORY
%
% I WROTE THIS PROGRAM TO HANDLE AN ANALYSIS NEED FOR A NOAA
% PROJECT APPLYING ITRDB TREE-RING CHRONOLOGIES TO CLIMATE
% STUDIES.  WE HAVE A FOXPRO DATABASE WITH ITRDB SITE AND 
% CHRONOLOGY INFORMATION, AND THE .CRN FILES AS IMPORTED OVER
% THE NET. THE OBJECTIVE IS TO BE ABLE TO READILY FORM
% NUMERIC MATRICES WITH CHRONS AS COLUMNS AND YEARS AS ROWS, 
% USING THE .CRN FILES FOR THE DATA, AND CONTROL FILES 
% FROM THE DATABASE AS SCREENING TOOLS.  
%
% THE APPROACH IS TO CULL CHRONOLOGIES FOR WHATEVER SUB-PERIOD
% AND OTHER CRITERIA, AND WRITE THE ANNUAL DATA FOR ALL 
% SELECTED CHRONOLOGIES AS A STRUNG-OUT VECTOR (SOV). A LINKING
% FILE WITH POINTERS TO ROWS WILL ALLOW LATER BUILDING OF 
% WORKABLE TIME SERIES MATRICES FROM THE SOV.   FILES 6 AND 7
% DESCRIBED ABOVE HOLD THE KEY OUTPUT.
%
% I originally did this work with a fortran program crn2sov.for
%
%
% BEWARE DIMENSIONS. CHECK FOR FOLLOWING VARIABLES. ASSUMING NUMBER
% OF SITES IS NSITES, AND MAXIMUM NUMBER OF YEARS IN ANY .CRN FILE IS
% MAXYRS
% 
%  YEARS(NSITES,2)
%  KYRS (NSITES),IC(NSITES),ID(NSITES),JSUM(NSITES)
%  X(MAXYRS), Y(MAXYRS)       


a = NaN;  % for convenience

% Input missing value code
title='Missing Value Code for Input Indices';
prompt = {'Enter the missing value code'};
def = {'9990'};
lineNo = 1;
msscode=inputdlg(prompt,title,lineNo,def);
msscode=str2num(msscode{1});  % as number

% Get name of file with input control
[file1,path1]=uigetfile('*.ctl','Input control file, file1');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

%--------  READ INFO OFF CONTROL FILE ----------------

% Get number of sites
c=fgetl(fid1);
nsites = str2num(c);

% Get specified range of years to be considered for including in sov
c=fgetl(fid1);
iyr1 = str2num(c(1:4));
iyr2=str2num(c(6:9));

% Name of file with list of names of input .crn files
c = fgetl(fid1);
pf2=strtok(c);

% Name of input file with list of start and end years of chronologies
c = fgetl(fid1);
pf3=strtok(c);

fclose(fid1);


%-------------  READ DATABASE LIST OF FIRST, LAST YEARS OF CHRONOLOGIES

fid3=fopen(pf3,'r');

% Want start year in col 1, end year in col 2
N = fscanf(fid3,'%5d',[2,nsites]);
N=N';

% Convert "8000" convention series to regular years
Ltemp=N(:,2)>2100;  % end year in database file greater than 2100AD
if any(Ltemp);
   N(Ltemp,1)=N(Ltemp,1)-8000;
   N(Ltemp,2)=N(Ltemp,2)-8000;
end
fclose(fid3);


%-------------- FROM DATABASE INFO, COMPUTE ROW INDEX OF CHRON IN SOV

% Offset of database start, end yrs of chrons from desired min start yr
% and max end yr of data in sov
idiff1=iyr1  - N(:,1);
idiff2=iyr2 - N(:,2);

% Revised start year -- start year of data to be used in sov.  Want
% to restrict data for sov to inclusive period iyr1 to iyr2
L1 = idiff1<0;
L2 = idiff2>0;

ic = repmat(iyr1,nsites,1);
id = repmat(iyr2,nsites,1);

if any(L1);  % chrons that start later than year iyr1
   ic(L1)=N(L1,1); 
end
if any(L2);  % chrons that end earlier than iyr2
   id(L2)=N(L2,2);
end

% Number of data values to be stored for each chron
nkeep = id-ic+1;

% Start row index in sov for each chron
ie = cumsum(nkeep)+1;
ie(nsites)=[];
ie = [1 ; ie];

% End row index
ig = ie + nkeep - 1;

% Store result
YRS = [(1:nsites)'  ic id ie ig];

% Initialize strungout vector
vlength=sum(nkeep);
v=a(ones(vlength,1),:);

%------ LOOP OVER THE CHRONOLOGIES, READING AND STORING DATA

% Open file of filenames
fid2 = fopen(pf2,'r');

% Loop
for n = 1:nsites;
   pfn = fgetl(fid2);
   pfn = strtok(pfn);
   igo = N(n,1); % start year of .crn file according to database
   isp = N(n,2); % end ...
   disp(n);
   disp(pfn);
   % get the vector of data and sample size for this seris
   [x,s,yr]=crn2vec2(pfn);
   
   % Check database start and end years against those in .crn file
   yron = yr(1);  yroff = yr(length(yr));
   if yron~=igo | yroff ~=isp;
      disp(['Site ' int2str(n)]);
      error('Database and .crn file start, end years inconsistent');
   end
   
   % Cull desired data segment for sov
   Lkeep = yr>=iyr1 & yr<=iyr2;
   xkeep = x(Lkeep);
   
   % Put the culled data in v
   v(ie(n):ig(n)) = xkeep;
   
end



function hcnpget1
% hcnpget1:  build .mat files of hcn monthly pcp and site info file
% CALL: hcnpget1.m
%
% Meko 10-2-97
%
%
%*********************  IN ******************************
%
% User is prompted for filenames of the station list 
% (typically stninf?.txt) and the main data file (hcnp.mat)
%     X in hcnp.mat holds the data for all 591 w. NA stations
%
%********************* OUT *****************************
%
% Individual .mat files are produced for each station in 
% stninf?.txt.  Files are put in a directory of choice, following
% a prompt  
%
% stnlst?.txt is same as stninf1.txt, but with sequential number in
% col 1
%
%
%*************************  NOTES ********************
%
% See howto\hcncook.txt and hcnpget1.txt for description of the hcnp.mat
%
% The stninf?.txt file is generated in foxpro from hcnpcp.dbf
% using a query (find stninf?.qpr  for example)


% Allow for max of 1500 files
b5=blanks(8);
CN = repmat(b5,1500,1); % will hold names of .mat files 

b6=blanks(6);
ID = repmat(b6,1500,1); % will hold station ids

b94 = blanks(94);
TX = repmat(b94,1500,1); % will hold full station info


% Open and close a placebo file just to tell program where to store
% output .mat files
[file1,pathout]=uigetfile('*.*','Click any file, need dir for output .mat files');
clear file1;


% Open station list file
[file1,path1]=uigetfile('stninf*.txt','Station list file');
pf1=[path1 file1];
fid1=fopen(pf1,'r');



%************* Find out how many stations to get, store their ids and output file names

% Skip first 3 lines.  Should be 2 blank lines and a third line with
% col headers with 'SEQNO' in cols 3-7
c=fgetl(fid1); % first blank line
c=fgetl(fid1); % second blank line
c=fgetl(fid1);  % header line
if ~strcmp(c(3:8),'SEQ_NO');
	fclose all
	error('Third line of stninf file needs SEQ_NO in cols 3-8');
end



k1=1;
nfiles=0;
while k1;
	c=fgetl(fid1);
	if ~(feof(fid1) & length(c)<94);
      nfiles=nfiles+1;
      stnid=c(10:15); % station id
      stnid=strjust(stnid); % right justify
      % Change any leading blank to zero
      if isspace(stnid(1));
         stnid(1)='0';
      end
      ID(nfiles,:)=stnid;
      CN(nfiles,:)=['pc' stnid]; % .mat file prefix
      c(10:15)=stnid;
		TX(nfiles,:)=c(1:94); % full line of site info
   else
      k1=0;
	
	end
end
CN=CN(1:nfiles,:);
ID=ID(1:nfiles,:);
TX=TX(1:nfiles,:);
yrsx = zeros(nfiles,2);  % will later hold start, end year for each station 
rows = zeros(nfiles,2); % will later hold start, stop row of data in X 
years = zeros(nfiles,2) ; % will hold database info on what should be the years

fclose(fid1);



%******************** BUILD MATRIX OF STARTING AND ENDING ROWS OF SERIES IN X
clc
disp('BUILDING MATRIX OF START AND END ROWS, yrs, IN X FOR DESIRED SERIES');
for n = 1:nfiles;
   c = TX(n,:);
   gorow = str2num(c(84:88));
   sprow = str2num(c(90:94));
   years(n,1)=str2num(c(73:76));
   years(n,2)=str2num(c(78:81));
   rows(n,:) = [gorow sprow];
end



%***********************  BUILD THE .MAT FILES FOR EACH STATION 
disp('STARTING TO MAKE .MAT STATION FILES');

% Load the mother monthly pcp file
[file2,path2]=uigetfile('hcnp.mat','Input file with all w N. Am data');
pf2=[path2 file2];
eval(['load ' pf2]); % data is in X

% Loop over stations
csave=blanks(6);
for n = 1:nfiles;
    id = ID(n,:);
    name=CN(n,:);
	 name2 = TX(n,17:26);
    disp(['   Getting data for: ' name ': ' name2]);
    rowgo=rows(n,1); % start row of station data in X
    rowsp=rows(n,2); % end row of station data in X
    x = X(rowgo:rowsp,:);  % station id, year, and 12 values
    
    % Check that id in col 1 is identical for all rows for this series
    idx = x(:,1);
    idstart = idx(1);
    if ~all(idx == idstart);
       disp(['Problem with station: ' name ' '  name2]);
       error('ids in first col of data not the same for all years');
    end
    
    % Check that the id from the data file matches the id from the database
    % First might need to slap a leading zero
    if idstart<100000;
       idtxt = ['0' int2str(idstart)];
    else
       idtxt = int2str(idstart);
    end
    
    if ~strcmp(idtxt,id);
       disp(['Problem with station: ' name ' '  name2]);
       disp(['id from hcnp.dbf:  ' id]);
       disp(['id from col 1 of data in X:  ' idtxt]);
       error('id from database does not match id from X (from hcnp.mat)');
    end
    
    % Drop the leading column (id) from x
    x(:,1)=[];
    
    % Check that the start and end year of x match that what database says
    yrx = x(:,1);
    yearson = years(n,1);
    yearsof = years(n,2);
    if min(yrx)~=yearson | max(yrx) ~= yearsof;
       disp(['Problem with station: ' name ' '  name2]);
       error('Years for this station not same as in database');
    end
    
    % Make sure years continuous
    d = diff(yrx);
    if ~all(d==1);
       disp(['Problem with station: ' name ' '  name2]);
       error('Year column not incrementing by one');
    end
    

    % Build filename for mat file
    pf3=[path2 name];
    
    % Rename x as Z and save in .mat file
    Z=x;
    eval(['save ' pf3 ' Z;']);
    
         
end



%**********  Build augmented station info

[file5,path5]=uiputfile('stnlst*.txt','Augmented station info file');
pf5=[path2 file5]; % will want this with the mat files
fid5=fopen(pf5,'w');


for n =1:nfiles
	str1=sprintf('%4.0f',n);
	str2=TX(n,:);
   str4 = [str1 str2];
   str4=str4(1:85);
	fprintf(fid5,'%s\n',str4);
end
fclose(fid5);


%**********  BUILD ASCII FILE OF STATION NAMES -- USEFUL IF LATER RUNNING 
% REGCLI2.M;  and ascii file of long-lat

[file6,path6]=uiputfile('filelst*.txt','List of .mat files for regcli2.m');
pf6=[path2 file6]; % will want this with the mat files
fid6=fopen(pf6,'w');


[file7,path7]=uiputfile('lonlat*.dat','Lon-lat file for regcli2.m');
pf7=[path2 file7]; % will want this with the mat files
fid7=fopen(pf7,'w');

for n =1:nfiles
   siteno=sprintf('%4.0f',n);
   fprintf(fid6,'%s  %s\n',CN(n,:),siteno);
   lon1 = str2num(TX(n,60:66));
   lat1 = str2num(TX(n,54:58));
   fprintf(fid7,'%7.2f %5.2f\n',lon1,lat1);
end

fclose all




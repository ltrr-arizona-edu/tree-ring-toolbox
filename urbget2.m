function urbget2
% urbget2:  get .mat files of urban heat island adjusted hcn monthly tmp
% CALL: urbget2
%
% Meko 10-1-97
%
%
%*********************  IN ******************************
%
% User is prompted for filenames of the station list 
% (typically stninf?.txt) and the main data file (urban.mea)
%
%********************* OUT *****************************
%
% Individual .mat files are produced for each station in 
% stninf?.txt.  Files are put in a directory of choice, following
% a prompt  
%
% stnlst?.txt is same as stninf1.txt, but with sequential number and
% with start and end year of data
%
%
%*************************  NOTES ********************
%
% See howto\hcntget1.txt for description of the urban.mea file
%
% The stninf?.txt file is generated in foxpro from hcntmpu.dbf
% using a query (find stninf?.qpr  for example)


% Allow for max of 1500 files
b5=blanks(8);
CN = repmat(b5,1500,1); % will hold names of .mat files 

b6=blanks(6);
ID = repmat(b6,1500,1); % will hold station ids

b71 = blanks(71);
TX = repmat(b71,1500,1); % will hold full station info


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
if ~strcmp(c(3:7),'SEQNO');
	fclose all
	error('Third line of stninf file needs SEQNO in cols 3-7');
end



k1=1;
nfiles=0;
while k1;
	c=fgetl(fid1);
	if ~(feof(fid1) & length(c)<71);
      nfiles=nfiles+1;
      ID(nfiles,:)=c(9:14); % station id
      CN(nfiles,:)=['tu' c(9:14)]; % .mat file prefix
		TX(nfiles,:)=c(1:71); % full line of site info
   else
      k1=0;
	
	end
end
CN=CN(1:nfiles,:);
ID=ID(1:nfiles,:);
TX=TX(1:nfiles,:);
yrsx = zeros(nfiles,2);  % will later hold start, end year for each station 

fclose(fid1);

%***********************  MAKE .DAT ASCII FILES 
disp('STARTING TO MAKE .DAT ASCII FILES');

% Open data file
[file2,path2]=uigetfile('urban.mea','Input data file');
pf2=[path2 file2];
fid2=fopen(pf2,'r');

csave=blanks(6);
k1=1;
kprev=0; % Last line read was not a match
for n = 1:nfiles;
    id = ID(n,:);
    name=CN(n,:);
	  name2 = TX(n,16:25);
    disp(['   Getting data for: ' name ': ' name2]);
   pf3=[path2 name '.dat'];
   fid3=fopen(pf3,'w');
      if n>1;
      if strcmp(csave(1:6),id);
         fprintf(fid3,'%s\n',csave(8:95));
      end
   end
   
   k1=1;
   while k1;
      c=fgetl(fid2);
      if feof(fid2);
         fclose(all);
         k1=0;
         break
      end
      
      cid = c(1:6);
      if strcmp(cid,id); % a match
         fprintf(fid3,'%s\n',c(8:95));
         kprev=1;
      else; % not a match
         csave=c; % in case this line is match for next file
         if kprev==1;  % last line was a match; so have finished a file
            fclose(fid3);
            kprev=0;  
            k1=0;
         else; %  this line not a match, last line not a match
         end
      end
   end; % while k1
end; % for n=1:nfiles

fclose all;


%**************** LOAD .DAT ASCII FILES AND CONVERT TO .MAT FILES, WITH SOME
%  QUALITY CONTROL TO REPLACE -99.99 WITH NAN AND TO CHECK FOR CONTINUITY OF YEARS

disp('STARTING TO BUILD .MAT FILES');

anan=NaN;
for n=1:nfiles;
   name = CN(n,:);
   disp(['   Working on : ' name]);
   if exist(name)==1;
     error('Bad coincidence: file name is also a variable name in this function');
   end
   
   pf3=[path2 name '.dat'];
   eval(['load ' pf3 ';']);
   eval(['Z = ' name ';']);
   
   % Replace -99.99 with NaN
   yr = Z(:,1);
	yrsx(n,:)=[min(yr) max(yr)];
   L=Z==-99.99;
   if any(any(L))
      Z(L)=anan;
   end
   
   % Make sure years monotonic increasing
   d = diff(yr);
   if ~all(d>=1);
      error('Years not montonic increasing');
   end
   
      
   % Insert any totally missing years with all-NaN data
   if ~all(d==1); % year does not increment by 1
      % Build a 13-col NaN matrix with year as col 1
      yrgo=min(yr); yrsp=max(yr);
      ntot = yrsp-yrgo+1;
      yrnew = (yrgo:yrsp)';
      Znew = repmat(NaN,ntot,13);
      Znew(:,1) = yrnew;
      % Make an index to rows of the new matrix telling where to insert data from Z
      irow = yr - yrgo+1;
      % Insert the data
      Znew (irow,2:13) = Z(:,2:13);
      Z=Znew;
   end
   
   
   
   
   
   % Save file
   eval(['save ' path2 name ' Z']);
  
   % Delete the ascii file
   eval(['delete ' pf3]);
   
end

%**********  Build augmented station info

[file5,path5]=uiputfile('stnlst*.txt','Augmented station info file');
pf5=[path2 file5]; % will want this with the mat files
fid5=fopen(pf5,'w');


for n =1:nfiles
	str1=sprintf('%4.0f',n);
	str2=TX(n,:);
	str3 = sprintf('  %4.0f %4.0f',yrsx(n,1),yrsx(n,2));
	str4 = [str1 str2 str3];
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
   siteno = sprintf('%4.0f',n);
   fprintf(fid6,'%s  %s\n',CN(n,:),siteno);
   lon1 = str2num(TX(n,59:65));
   lat1 = str2num(TX(n,53:57));
   fprintf(fid7,'%7.2f %5.2f\n',lon1,lat1);
end




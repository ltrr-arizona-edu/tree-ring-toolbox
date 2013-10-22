function urbanget
% urbanget:  get .mat files of urban heat island adjusted hcn monthly tmp
% CALL: urbanget
%
% Meko 6-11-97
%
%
%*********************  IN ******************************
%
% User is prompted for filenames of the station list 
% (typically stnlist.txt) and the main data file (urban.mea)
%
%********************* OUT *****************************
%
% Individual .mat files are produced for each station in 
% stnlist.txt.  
%
%
%*************************  NOTES ********************
%
% See howto\hcntget1.txt for description of the stnlist.txt and
% urban.mea files.


% Allow for max of 1500 files
b5=blanks(5);
CN = repmat(b5,1500,1); % will hold names
b6=blanks(6);
ID = repmat(b6,1500,1); % will hold station ids



% Open station list file
[file1,path1]=uigetfile('stnl*.txt','Station list file');
pf1=[path1 file1];
fid1=fopen(pf1,'r');



%************* Find out how many stations to get, store their ids and output file names

k1=1;
nfiles=0;
while k1;
	c=fgetl(fid1);
	if ~feof(fid1) & length(c)>50;
      nfiles=nfiles+1;
      ID(nfiles,:)=c(1:6);
      i1 = find(isspace(c));
      i1=i1(1);
      c(1:i1)=[];
      c=fliplr(c);
      c=deblank(c)
      c=fliplr(c);
      CN(nfiles,:)=c(1:5); % store first part of output file name
   else
      k1=0;
	
	end
end
CN=CN(1:nfiles,:);
ID=ID(1:nfiles,:);

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
    disp(['   Getting data for: ' name]);
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

disp('STARTING TO BUILE .MAT FILES');

anan=NaN;
for n=1:nfiles;
   name = CN(n,:);
   disp(['   Working on : ' name]);
   if exist(name)==1;
     error('Bad coincidence: file name is also a variable name in this function');
   end
   
   pf3=[path2 name '.dat'];
   eval(['load ' pf3]);
   eval(['Z = ' name]);
   
   % Replace -99.99 with NaN
   yr = Z(:,1);
   L=Z==-99.99;
   if any(any(L))
      Z(L)=anan;
   end
   
   % Make sure years continuous
   d = diff(yr);
   if ~all(d==1);
      error('Year column not incrementing by one');
   end
   
   % Save file
   eval(['save ' path2 name ' Z']);
  
   % Delete the ascii file
   eval(['delete ' pf3]);
   
end

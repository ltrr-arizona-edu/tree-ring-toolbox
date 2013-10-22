function rwfiles
% rwfiles: cell matrix of prefixes of rw (or lww or eww( filenames in current directory
% CALL: rwfiles;
%
% Meko 3-17-98; mod 2-8-99
% 
%********************* OUT ********************
%
% No output args.
%
% Builds cell matrix of filename prefixes and stores as specified 
% variable in a specified .mat file.  Names of variable and filename
% prompted for.
%
% Optinally, a corresponding ascii file, put in same directory
%   as the input files

% Prompt for name of cell variable to store prefixes
tit='Name of cell variable to store the prefixes';
prompt={'Enter name:'};
def={'rwpre'};
lineNo=1;
answer=inputdlg(prompt,tit,lineNo,def);
cellnm=answer{1};

% Prompt for name of .mat file to store cell variable with prefixes
tit='Name of .mat file to store cell vbl of prefixes';
prompt={'Enter name of .mat file:'};
def={'prefixes'};
lineNo=1;
answer=inputdlg(prompt,tit,lineNo,def);
fn=answer{1};

% Prompt for file suffix of desired files 
prompt={'Enter suffix of files to be included in list:'};
def={'.rw'};
tit1='Suffix for Desired Files';
lineNo=1;
answer=inputdlg(prompt,tit1,lineNo,def);
suff=answer{1};
lensuff = length(suff); % length of suffix, including the .

% Initialize the cell variable as empty
eval([cellnm ' = [];']);

% Get names of files in current directory and store structure of info in D
D=dir;  % puts file info into structure D
nfls1 = size(D,1);  % number of files, including directories, in current dir
C=struct2cell(D); % convert structure to cell. C(1,:) holds file names
c = C(1,:); % file names, in cell array

c1 = sortrows(char(c)); % convert cell array c to char array, sorted alphabetically

ncount=0; 
% Loop over filenames
for n = 1:nfls1;
   d1 = c1(n,:); %  a single filename, as char -- could also be a directory
   d1 = strtok(d1); % truncate trailing blanks
   i1 = findstr(d1,lower(suff));
   i2 = findstr(d1,upper(suff));
   if ~(isempty(i1) & isempty(i2));
      if isempty(i1);
         i1=i2;
      else
         i1=i1;
      end
      i1=i1(length(i1)); % to make sure pointing to start of last occurrence of suff
      % Make sure suffix contains exactly the right number of characters
      if length(d1)==i1+lensuff-1;
         ncount=ncount+1;
         d1 = strtok(d1,'.'); % drop off the suffix
         eval([cellnm '{ncount}=d1;']);
         
      else
      end
   end
   
end

% Error if no files with the desired suffix
eval(['Lnone = isempty(' cellnm ');']);
if Lnone;
   error(['No files found with suffix: ' suff]);
end

% Make output file with cell variable of file names
eval(['save ' fn ' ' cellnm ';']);




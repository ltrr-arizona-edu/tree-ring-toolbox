function pre=dirfls3
% dirfls3: cell matrix of prefixes of filenames in current directory
% CALL: pre=dirfls3;
%
% Meko 3-17-98
% Builds cell matrix of prefixes of filenames in  current directory
% 
%********************* OUT ********************
%
% pre:  cell matrix of filename prefixes
% Optinally, a corresponding ascii file, put in same directory
%   as the input files



pre=[];

D=dir;  % puts file info into structure D
nfls1 = size(D,1);  % number of files, including directories, in current dir
C=struct2cell(D); % convert structure to cell. C(1,:) holds file names
c = C(1,:); % file names, in cell array

prompt={'Enter suffix:'};
def={'.rw'};
tit1='Suffix for Desired Files';
lineNo=1;
answer=inputdlg(prompt,tit1,lineNo,def);
suff=answer{1};
lensuff = length(suff); % length of suffix, including the .

c1 = char(c); % convert cell array c to char array

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
         pre{ncount}=d1;
         
      else
      end
   end
   
end

if isempty(pre);
   error(['No files found with suffix: ' suff]);
end



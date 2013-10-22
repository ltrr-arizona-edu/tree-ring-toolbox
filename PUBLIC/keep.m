function keep(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
% KEEP is complementary to the "clear" command:
% instead of clearing the supplied variables, it
% keeps them and deletes the rest. 
%
% For example, the command  " keep x y " clears
% all variables from the workspace except x and y.

% Written by Eric Breitenberger. Please send any bug
% reports or other comments to: eric@gi.alaska.edu.

oldfile=get(0,'DiaryFile'); % get current DiaryFile
diaryfile='keep.txt';
% If the diary file already exists, delete it
if (eval(['exist(''' diaryfile ''')']))==2
  delete(diaryfile)
end

% Jump up to main workspace, write variable names to diary file,
% and jump back down. Note the use of 'eval' to pass the diaryfile 
% argument back out to the main workspace.
eval(['dbup, diary(''' diaryfile ''')'])
who
diary off
dbdown

% Read the diary file, then delete it
fid=fopen(diaryfile,'r');
filestr=fread(fid)';
fclose(fid);
delete(diaryfile)
set(0,'DiaryFile',oldfile) % restore old DiaryFile

% Find the colon (58) in "Your variables are:" and truncate to that point:
ind=find(filestr==58);
filestr=filestr(ind(1)+1:length(filestr));

% Find all carriage returns and line feeds and convert to spaces (32).
% PCs use CR-LF (13 10), Unix boxes use LF (10), Macs use CR (13).
ind=find(filestr==13); % Carriage returns
ind=[ind find(filestr==10)]; % Line feeds
filestr(ind)=32*ones(1,length(ind));

% Now go through and change multiple spaces to single spaces:
filestr=[32 filestr 32]; % make sure to start and end with a space                  
index=find(filestr==32);
spaces=zeros(1,length(filestr));
spaces(index)=ones(1,length(index)); % 'spaces' =1 for spaces, =0 for any other char.
d=diff([1 spaces]);
endvar=find(d==1); % these are the indices of the first space after each var
if endvar==[], error('No variables in workspace!'), end
filestr(endvar)=999*ones(1,length(endvar)); % Change the first space after each var to 999
index=find(filestr==32); 
filestr(index)=[]; % strip spaces
endvar=find(filestr==999); % add a single space at beginning
numvars=length(endvar);   % number of workspace variables
filestr(endvar)=32*ones(1,numvars); % put single spaces back in for the 999s
% Filestr now contains all the workspace variables, separated by trailing spaces
filestr=[32 filestr]; % need that initial space again

% For each input argument, find if it exists in filestr, and where it sits:
remove=[]; % this will hold the indices to be removed from filestr
numremove=0; % number of removed variables
notfound=0;  % number of inputs not found in filestr
for i=1:nargin
  temp=eval(['a' num2str(i)]);
  var=[32 abs(temp) 32];   % space pad to ensure unique strings
  k=findstr(var,filestr);
  if isempty(k)
    notfound=notfound+1;
    disp(['Warning: Input argument not found in workspace:' var])
  else 
    remove=[remove k+1:k+length(var)-1]; 
    numremove=numremove+1;
  end 
end

filestr(remove)=[]; % remove the found variables from filestr

% set up string(s) to pass to 'clear': If there are more than 50 variables,
% the 'clear' command crashes, so variables will be passed 50 at a time.
numclear=numvars-numremove;
if numclear<=50
  eval(['dbup, clear ' filestr ', dbdown'])  % Clear the unwanted stuff:
else % Use the same code as up above to find the variable locations
  index=find(filestr==32);
  spaces=zeros(1,length(filestr));
  spaces(index)=ones(1,length(index));
  d=diff([1 spaces]);
  stvar=find(d==-1); % these are the indices of the first char. in each var
  endvar=find(d==1); % these are the indices of the first space after each var
  numcalls=floor((numclear-1)/50)+1;  % compute number of calls to 'clear'
  for k=1:numcalls % for each call, get the start and endpoints, then clear
    istart=50*(k-1)+1;
    iend=min(numclear,50*k);
    string=filestr(stvar(istart):endvar(iend)-1);
    eval(['dbup, clear ' string ', dbdown'])
  end
end


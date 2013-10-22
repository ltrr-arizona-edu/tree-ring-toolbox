function sov2mon1
% sov2mon1:  steve leavitt single-column ascii monthly pcp or tmp to 13-col Z in .mat file
% CALL: sov2mon1
%
% Meko 10-3-98
%
%************* IN 
%
% No args
%
% User prompted to pick an input  file, which is either .pcp or .tmp
% User prompted for start and end years of the data
%
%*********** OUT
%
% Data automatically output  as 13-col Z in .mat file with same prefix as input file,
% in same directory
%
%********** NOTES
%
% Pcp assumed in inches.  Tmp in deg F
% function checks that number of values is 12 times number of expected years
% User ideally changes to the directory that the input files are in before
%       running this function

ktype = menu('Data of Which Type?','pcp','tmp');
if ktype==1;
   flspec='*.pcp';
else;
   flspec='*.tmp';
end

[file1,path1]=uigetfile(flspec,'Input monthly climate data file');
pf1=[path1 file1];


% Load ascii input file
eval(['load ' pf1]);
name1 = strtok(file1,'.');
eval(['X = ' name1 ';']);

% Check that suffix matches expected data type
str1=flspec(3:5);
len1 = length(pf1);
str2 = lower(pf1((len1-2):len1));
if ~strcmp(str1,str2);
   error(['Wanted ' flspec ';  got ' pf1]);
end

titdlg='Expected Years of Input Monthly Data';
prompt = {'Start year','End year'};
def = {'1980','1995'};
lineNo=1;
answer=inputdlg(prompt,titdlg,lineNo,def);

yrgo = str2num(answer{1});
yrsp = str2num(answer{2});
nyr = yrsp-yrgo+1;
nrow1 = nyr*12; % expected number of rows of strung out data

% Check for correct number of values
[mX,nX]=size(X);
if nX~=1;
   error('X is not a cv');
end
if mX ~= nrow1;
   disp(pf1);
   disp(['nrow1 = ' int2str(nrow1)]);
   disp(['N rows in X = ' int2str(mX)]);
   error('nrow1 should equal mX (see above)');
end

%----- RESHAPE

Y = (reshape(X,12,nyr))';
yr = (yrgo:yrsp)';
Z = [yr Y];  % add year column

%----- SAVEW

name2 = lower(name1);
eval(['save ' name2 ' Z;']);


pf1;
 


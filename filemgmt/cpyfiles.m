function cpyfiles
% cpyfiles: copy files from various directories to a specified directory
% CALL: cpyfiles
%
% Meko 4-8-99
%
% Why?  To copy m-files for class use to a directory for testing, or to 
% a destination such as a zip drive or network directory
%
%********** IN ****
%
% No input args
% User prompted for information:
%   -file with source drive/path/name of files desired
%   -via dialog for drive/path/ of destination
%
%******* OUT 
%
% No output args
% Action is that the files have been copied to the destination
% Format of input file: one drive\path\filename per line


kmen1 = menu('Choose default situation',...
    'Testing geos595e scripts',...
    'Testing treetool',...
    'Testing rwmeas');
if kmen1==1;
    defsrc = 'c:\geos595e\mfiles.txt';
    defdest= 'c:\workclass\';
elseif kmen1==2;
    defsrc= 'c:\mlb\filemgmt\flstree.txt';
    defdest = 'c:\mlb\treetool\';
elseif kmen1==3;
    defsrc= 'c:\mlb\velmex\flvelmex.txt';
    defdest = 'c:\work6\';
    
end;

% Prompt for name of file with list of filenames
Prompt='Input Filename:';
Title='Input file with list of files to be copied';
LineNo=1;
DefAns={defsrc};
Answer = inputdlg(Prompt,Title,LineNo,DefAns);
A=char(Answer{1});

% Prompt for destination drive\path\ for files
Prompt='Target drive\directory:';
Title='Drive and directory to which files are to be copied';
LineNo=1;
DefAns={defdest};
Answer = inputdlg(Prompt,Title,LineNo,DefAns);
B=char(Answer{1});


% Read names of files to be copiedinto a column-cell of strings
F = textread(A,'%s','delimiter','\n','whitespace','');

% Calc number of files
nfile = size(F,1);

clc;
disp([int2str(nfile) ' files to be copied']);

% Loop over the files to be copied
for n = 1:nfile
   f=F{n}; % a filename, which may or may not contain drive\directory
   
   disp(f);
   
   % Cull the filename part off of drive\directory
   islash = findstr(f,'\');
   if isempty(islash);
      fname=f;
   else
      fname=f;
      fname(1:(max(islash)))=[];
   end
   
   % Prepend the destination drive\directory
   fname = [B fname];
   
   % Copy the file from source to destination
   [status,msg]=copyfile(f,fname);
   if status==0;
      disp(msg);
      error(' ');
   end
   
   
end
msg1=[int2str(nfile) ' files successfully copied'];
clc;
disp(msg1);



function treecpy2
% treecpy1: copy crvfit and related standardization functions to a single directory
% treecpy2;
% Last revised 9-9-99
%
% treecpy2.m is needed to gather tree-ring toolbox functions from various 
% subdirectories under \mlb\ into a single directory. Used to facilitate 
% gathering the functions for testing on Kaimeis PC
%
%*** INPUT -- no arguments
%
% User prompted to point to a .txt file holding path\name of functions
% User prompted to name directory to copy functions to
%
%*** OUTPUT --- no arguments
%
% Action is to copy the files to the specified target directory
%
%
%*** REFERENCES
%*** UW FUNCTIONS CALLED
%*** NOTES

%--- Get .txt files with function names
[file1,path1]=uigetfile('c:\manip\stzfls.txt','Input file of names of functions to copy');
pf1=[path1 file1];
flnms = textread(pf1,'%s','delimiter','\n','whitespace',''); % cell of function names

%-- Get target directory
prompt={'Enter the target directory:'};
def={'d:\kaimei'};
dlgTitle='Target directory to copy toolbox functions to';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
path2=answer{1};

%-- Compute number of functions to copy
nfuns = size(flnms,1);

%-- Loop over functions
disp('Files copied:');
for n = 1:nfuns;
   fl1 = flnms{n};
   % Pull the filename off the path
   i1=findstr(fl1,'\');
   if isempty(i1); % if no slashes
      % no action needed; fl1 is the file name
      fl2=fl1;
   else;
      fl2=fl1;
      fl2(1:max(i1))=[];
   end
   
   % Build the target path\filename
   pf2 = [path2 '\' fl2];
   
   % Copy the file
   copyfile(fl1,pf2);
   
   disp([fl1 ' to ' pf2]);
   
   
end

   
      






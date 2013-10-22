function copyrw1
% copyrw1:  copy specified .rw files from one directory to another
% copyrw1;
% Last revised 2-22-00
%
% Used to copy selected rw files to a target directory
%
%*** INPUT -- NO ARGS
%
% You have built an ascii .txt file (e.g., jac98.txt) with a .rw filname on
% each row.   Example:
% jac45a.rw
% jac03b.rw
%
% You are prompted for the name of this file (*.txt) 
%
% You are also prompted for the path\filename of the target directory
%
%*** OUTPUT -- NO ARGS
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%*** NOTES
%
% The working directory should be the source directory of the .rw files

% PROMPT FOR MODE OF COPY
kmen1=menu('Choose one',...
    'Names in filename list have no . or suffix; all are are .rw',...
    'Names in filename list have . and suffix (e.g. fry05.LWW)');


% PROMPT FOR INPUT-CONTROL FILE
[file1,path1]=uigetfile('*.txt','File with list of .rw files to be copied');
pf1=[path1 file1];
file = textread(pf1,'%s','delimiter','\n','whitespace','');
nfiles=size(file,1);

% PROMPT FOR TARGET DIRECTORY 
prompt={'Enter the target directory:'};
def={'c:\work\'};
dlgTitle='Directory .rw files are to be copied to';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
path2=answer{1};


% COPY THE FILES
for n = 1:nfiles;
    if kmen1==1;
        src=[file{n} '.rw'];
    elseif kmen1==2;
        src=[file{n}];
    end;
   dest=[path2 src]; 
   [status,msg] = copyfile(src,dest);
   if status==0;
      error(msg);
   else;
      disp([src ' sucessfully copied']);
      
   end;
end;









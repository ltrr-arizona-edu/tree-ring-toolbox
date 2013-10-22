function rwlcull
% rwlcull:  build a .rwl from ring width series culled from a source rwl
% rwlcull;
% Last revised 8-14-00
%
% Builds a .rwl file from ring widths culled from a source .rwl file. WritteN
% to edit .rwl file so that only series used in chronology are included in
% the .rwl sent to ITRDB
%
%*** INPUT -- NO ARGS
% 
% Your are prompted for the name of an ascii .txt 
% file (e.g., jaclist2.txt) with the list of
% desired ring-width series.   Example:
% jac03b
% jac43x
%
% You are prompted for the name of the source .rwl file
% to cull the ring-width series from (e.g., 'nv060.rwl') 
%
% You are prompted for the path and filename of the .rwl file
% to contain the culled series
%
%*** OUTPUT -- NO ARGS
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% The series names in the .rwl file are assumed to be all capital letters, 
% but you may make the list with small letters because the list is made 
% upper case before the comparison.
%
% The working directory should be the directory with the source .rwl file

% PROMPT FOR LIST WITH DESIRED RING-WIDTH SERIES
[file1,path1]=uigetfile('*.txt','Infile with list of desired series for the target .rwl file');
pf1=[path1 file1];
file1 = textread(pf1,'%s','delimiter','\n','whitespace','');

nfiles=size(file1,1);
k=file1{nfiles};
if all(isspace(k));
    file1(nfiles)=[];
    nfiles=nfiles-1;
end;

% PROMPT FOR SOURCE .RWL FILE
[file2,path2]=uigetfile('*.rwl','Source .rwl file');
pf2=[path2 file2];
file2 = textread(pf2,'%s','delimiter','\n','whitespace','');
nlines=size(file2,1);

% PROMPT FOR TARGET DIRECTORY 
prompt={'Target directory:','Target .rwl file'};
def={'c:\work\','joker.rwl'};
dlgTitle='Directory/file of targert .rwl file';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
path3=answer{1}
file3=answer{2};
pf3=[path3 file3];


% LOOP OVER THE FILENAMES OF THE DESIRED RING WIDTHS
L = logical(zeros(nlines,1));

for n = 1:nfiles;
   rwname=file1{n};
   rwname=upper(rwname);
   nchar=length(rwname);
   L1=strncmp(file2,rwname,nchar); % mark rows with desired series
   n1 = sum(L1);
   if n1==0;
      error([rwname ' not in ' pf2]);
   else;
      L(L1)=1;
      disp([rwname  ' successfully copied']);
   end;
end;

A=char(file2(L));
[mA,nA]=size(A);
fid1=fopen(pf3,'w');
for n = 1:mA;
   a=A(n,:);
   fprintf(fid1,'%s\n',a);
end;
fclose(fid1);
 







disp('here');











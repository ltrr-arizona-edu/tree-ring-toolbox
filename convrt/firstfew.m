function firstfew
% firstfew:  write first few lines of a file to another file
% firstfew;
% Last revised 4-19-00
%
% Needed because some alphanumeric data files are too huge for testing program code,
% and wordpad, edix, etc warp lines with line wrap
%


[file1,path1]=uigetfile('*.dat','Input big file');
pf1=[path1 file1];
[file2,path2]=uiputfile('*.*','Output subset of big file');
pf2=[path2 file2];

fid1=fopen(pf1,'r');
fid2=fopen(pf2,'w');

prompt={'Enter number of lines:'};
def={'500'};
dlgTitle='Number of lines to transfer';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
nkeep=str2num(answer{1});

for n=1:nkeep;
   c=fgetl(fid1);
   fprintf(fid2,'%s\n',c);
end;
fclose(fid1);
fclose(fid2);
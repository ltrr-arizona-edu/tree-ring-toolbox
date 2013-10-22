function S=dirfls(dir_name,suffix,kopt)
% dirfls: string matrix and ascii list of files in a directory
% CALL: S=dirfls(dir_name,suffix,kopt);
%
% Meko 3-17-97
% Builds ascii file of names of files in directory dir_name.
% If two input args, considers only files with suffix "suffix"
%
%**************** IN ***************************
%
% dir_name (1 x ?)s directory name
% suffix (1 x ?)s  suffix (without the period)
% kopt(1 x 1)i 
%   ==1 make ascii file
%   ==2 do not make ascii file
% Example:  dirfls('c\jack',suff,1), where suff is 'rwl'
%
%********************* OUT ********************
%
% S string matrix of file names
% Optinally, a corresponding ascii file, put in same directory
%   as the input files

kmen1=menu('Choose case of filenames',...
    'Upper',...
    'Lower');
if kmen1==1;
    casefn='Upper';
else;
    casefn='Lower';
end;



D=dir(dir_name);

[m1,n1]=size(D); % m1 is number of files in directory

S=blanks(8);
str2=['.' suffix];
% Loop over files
for n = 1:m1;
   dname=getfield(D(n),'name');
   % findstr tries to find occurrences of shorter string in longer.  Must make sure
   % str2 is not the longer of the two strings
   len1=length(dname);
   len2=length(str2);
   if len1<=len2;
      nblank=len2-len1+1;
      dname=[dname blanks(nblank)];
   end
   
   if strcmp(casefn,'Upper');
        k=findstr(upper(dname),str2);
    else;
        k=findstr(lower(dname),str2);
    end;
   if ~isempty(k);
      S=str2mat(S,dname);
   
   end
   
end
S(1,:)=[];

% Make ascii file
[ms,ns]=size(S);
if kopt==1;
    [file1,path1]=uiputfile('*.txt','Output list of file names');
    pf1=[path1 file1];
    fid1=fopen(pf1,'w');
    fmt1='%s\n';
    for k=1:ms;
        s=deblank(S(k,:));
        fprintf(fid1,fmt1,s');
    end
elseif kopt==2;
else
    error('kopt must be 1 or 2');
end
 
fclose all

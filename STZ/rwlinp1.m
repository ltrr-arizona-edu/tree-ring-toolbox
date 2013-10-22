function rwlinp1
% rwlinp1: .rwl files converted to sov .mat files of ring-width series
% rwlinp1.m is basically a driver function that uses rwlinp.m to convert lots of
% .rwl files at one time.  The user can optionally specify separate directories
% for finding the .rwl files and putting the .mat files and .tmp files.  Or, in
% a simple case, the user might just have all kinds of files in one directory.
%
% Meko 3-17-97, to convert san pedro river basin tree-ring  data .rwl files 
%
% 
%
%
% Open file with list of .rwl files
[file1,path1]=uigetfile('*.txt','File with names of .rwl files');
pf1=[path1 file1];
fid1=fopen(pf1,'r');


% Initial read of fid1 to count number of files
nfiles=0;
k1=1;
while k1;
   c=fgetl(fid1);
   if feof(fid1);
      k1=0;
   else
      nfiles=nfiles+1;
   end
end; % of while k1


% Rewind the file list
frewind(fid1);


% Ask whether want to specify .rwl directory rather than read off .rwl filenames
ButtonName=questdlg('Do you want to separately specify the path?');
switch ButtonName,
case 'No'
   k2=1;
case 'Cancel'
case 'Yes'
   k2=2;
   prompt='Enter path, ending with backslash, no quotes';
   title='Path to Input .rwl Files';
   answer=inputdlg(prompt,title) ;
   pathin=answer{1};
   
   prompt='Enter path, ending with backslash, no quotes';
   title='Path to Output .mat and .tmp Files';
   answer=inputdlg(prompt,title) ;
   pathout=answer{1};
end


% Call rwlinp to convert files
for n = 1:nfiles
   disp(['File # ' int2str(n)]);
   flname=strtok(fgetl(fid1));
   if k2==1;
      rwlinp(flname);
   elseif k2==2
      rwlinp(flname,pathin,pathout);
   else
      error('Invalid k2 value');
   end
end

fclose all;


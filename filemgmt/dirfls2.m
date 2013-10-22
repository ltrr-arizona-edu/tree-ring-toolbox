function dirfls2
% dirfls2: modify file list of directory files by adding a path and changing suffix
%
% Meko 3-19-97
%

% Get file with list of files
[file1,path1]=uigetfile('*.txt','Input file list');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

% Get filename for output
[file2,path2]=uiputfile('*.txt','Output file list');
pf2=[path2 file2];
fid2=fopen(pf2,'w');


% Get the path to add to filenames
prompt='Enter the path, including the ending backslash';
title=' ';
lineNo=1;
def={'D:\jill\'};
pref=inputdlg(prompt,title,lineNo,def);
pref=char(pref);


action=1;



k1=1;
while k1;
   c=fgetl(fid1); % get a filename
   if ~feof(fid1);
      switch action
      case 1 ;% add path and remove suffix
		  % remove any trailing or leading blanks
		  c=deblank(c);
		  c=fliplr(c);
		  c=deblank(c);
		  c=fliplr(c);
         c=strtok(c,'.'); % remove suffix
         c=[pref c]; % add path
      case 2
         
      otherwise
      end
      fprintf(fid2,'%s\n',c);
   else
      k1=0
   end
end

fclose(fid1);
fclose(fid2);

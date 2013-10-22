function urbhist1
% urbhist1: cull one state's station history info from urban heat adjusted hcn station hist file
% CALL: urbhist1
%
% Meko 6-13-97
%
%***************  IN ****************************************
%
% User prompted for file name <station.his> previously downloaded from ncdc
% User prompted for 2-digit state code  for desired state
% 
%
%*************  OUT ************************
%
% User prompted for file name of culled station histories



%**********  Open
[file1,path1]=uigetfile('station.his','Input hcn US station history file');
pf1=[path1 file1];
fid1=fopen(pf1,'r');


prompt={'Enter two-digit state code'};
def={'02'};
title='Selection of State for Output History File';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);

state=answer{1};
txt1=['State ' state];

[file2,path2]=uiputfile('stnhist.txt',['Output station history file, ' txt1]);
pf2=[path2 file2];
fid2=fopen(pf2,'w');


k1=0;
% check first line in case want alabama
c=fgetl(fid1),
if strcmp(c(1:2),state);
   k1=1;
   fprintf(fid2,'%s\n',c);
end
k2=1;
while k2==1;
   c=fgetl(fid1);
   if strcmp(c(1:2),state);
      k1=1;
      fprintf(fid2,'%s\n',c);
   else;
      if k1==0; % have not yet hit the first line for the key state
         % no action needed
      else
         % have finished last line for key state
         fclose (fid2);
         k2=0;
      end
   end
end
fclose all;
   
      


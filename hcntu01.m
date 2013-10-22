function hcnt01
% hcnt01:  Make an ascii file of urban heat adjusted HCN tmp station info
% CALL: hcnt01
%
% Meko 10-1-97
%
%************** IN *******************************
%
% User prompted for names of two files:
%
% 1. input station inventory file as downloaded from NCDC (station.inv)
% 2. output ascii file of reformatted info suitable for excel-->foxpro
%
%*************** OUT *****************************
%
% The ascii output file of station info.  Column data delimited by blanks
%
%
%************** NOTES *************************
%
% 


% Open the input inventory file
[file1,path1]=uigetfile('*.inv','Station inventory, from NCDC');
pf1=[path1 file1];
fid1=fopen(pf1,'r');


% Open the output ascii file
[file2,path2]=uiputfile('*.txt','Output reformatted station inventory');
pf2=[path2 file2];
fid2=fopen(pf2,'w');


% Loop over the 1221  US HCN stations
for n = 1: 1221;
   c=fgetl(fid1);
   
   % sequential number of station in database
   str1=sprintf('%4.0f',n);
   
   % station id (6 character)
   str2 = sprintf('%s',(c(1:6)));
   
   % Name of mat file;  "tu" indicates, temperature, urban heat island adjusted
   str3 = ['tu' str2];
   
   % Climatic division
   str4= sprintf('%s',c(8:9));
   
   % Still active"
   ctemp=c(10);
   if strcmp(ctemp,'*');
      str5='F';
   else
      str5='T';
   end
   
   % State
   str6 = c(32:33);
   
   % Longitude
   ctemp=str2num(c(18:24));
   str7 = sprintf('%7.2f',ctemp);
   
   % Latitude
   ctemp=str2num(c(12:16));
   str8 = sprintf('%5.2f',ctemp);
   
   % Elevation (ft)
   str9 = c(26:30);
   
   % Station name
   str10 = c(35:64);
   
   stralla=[str1 ' ' str2 ' ' str3 ' '  str4 ' ' str5 ' '];
   strallb=[str6 ' ' str7 ' ' str8 ' ' str9 ' ' str10];
   strall=[stralla strallb];
   
   % Copy line to output file
   fprintf(fid2,'%s\n',strall);
   
end




fclose all;

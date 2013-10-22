function clicon6
% clicon6:  get division 0207 monthly data from downloaded ncdc div pdsi or related file
% CALL: clicon6
%
%




[file1,path1]=uigetfile('9705.*','Input data ftpd data file');
pf1=[path1 file1];

fid1=fopen(pf1,'r');


[file2,path2]=uiputfile('*.dat',['Outfile for ' file1]);
pf2=[path2 file2];
fid2=fopen(pf2,'w');

id='0207';


kon=0;
k1=1;
while k1
   
   if feof(fid1);
      k1=0;
      fclose all;
      break
   end
   
   c=fgetl(fid1);
   if length(c)<4;
      k1=0;
      break
   end
   
   %disp(c);
   c4=c(1:4);
   if strncmp(id,c4,4);
      kon=1;
      %disp(c);
      fprintf(fid2,'%s\n',c);
   else
      if kon==1;
         k1=0;
         fclose all;
         break
      end
      
   end
end


fclose all


function azcull1
% azcull: cull az lines from ncdc 19th century clim files
% CALL:  azcull1;
%
%**************** IN ***********************************8
%
% User is prompted for two file names. First is input file of 
% whatever downloaded over net.  Second is the file to append
% the az lines to.
%
% The input file must of course exist.  The output file
% will be created if it does not already exist.  Data in the
% input whose lines start with AZ will be appended to the output file


[file1,path1]=uigetfile('joke.dat','Input downloaded file');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

[file2,path2]=uiputfile('*.txt','Output file to append az lines to');
pf2=[path2 file2];
if exist(pf2)==2;
	fid2=fopen(pf2,'a');
else;
	fid2=fopen(pf2,'w');
end


k1=0; % flag that will be 1 when hit first AZ line
k2=1; % for while
while k2;
	c=fgetl(fid1);
	if feof(fid1);
		k2=0;
		break;
	else
		if strcmp(c(1:2),'AZ')
			k1=1;
			fprintf(fid2,'%s\n',c);
		else;
			if k1==1;
				k2=0;
				break;
			else;
				% still no first hit
			end
		end
	end
end

fclose all;
		




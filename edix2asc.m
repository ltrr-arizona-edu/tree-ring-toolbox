function edix2asc
%
% Word-wrapped edix file to ascii file. Undos problem encountered
% with hcn data.  I accessed an hcn divisional pcp data, more than
% 80 cols, and copied some stations' records to another edix
% window. But the data was wrapped in the new file, so that a year's
% data took 2 lines.  Not good.
%
% Meko 12-10-96


% get the file
[file1,path1]=uigetfile('*.dat',' Input wrapped file');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

% Open a write file
[file2,path2]=uiputfile('*.dat',' Output file');
pf2=[path2 file2];
fid2=fopen(pf2,'w');




k1=1;
while k1;
	c1=fgetl(fid1);
	if ~feof(fid1);
		c2=fgetl(fid1);
		c=[c1 c2];
		fprintf(fid2,'%s\n',c);
	else
		k1=0;
	end
end

fclose all

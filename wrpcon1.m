function wrpcon1(n)
%
% Convert word-wrapped HCN divisional pcp file by combining pairs of 
% lines.  Also cuts off the first n characters of each starting
% line, meaning gets rid of the id info before the year.
%
%
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
		c1(1:n)=[];
		c2=fgetl(fid1);
		c=[c1 c2];
		fprintf(fid2,'%s\n',c);
	else
		k1=0;
	end
end

fclose all

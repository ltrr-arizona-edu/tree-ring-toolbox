function cooplst1
%
% Scan coop station list as obtained from NCDC.  Write subset of rows
% for AZ and NM stations only to coop2.txt


file1=uigetfile('coop.txt','ncdc coop stn info file');
fid1=fopen(file1,'r');
file2=uiputfile('coop2.txt','state-culled stn info');
fid2=fopen(file2,'w');

% First skip initial lines.  Station info will start after 3rd 
% blank line
n1=0;
k1=1;
while k1==1;
	c=fgetl(fid1);
	if all(isspace(c));
		n1=n1+1
	end
	if n1==3; % have read 3 blank lines
		k1=0;
	end
end


k2=1
while k2==1;
	c=fgetl(fid1);
	%disp(['length(c) = ' int2str(length(c))  ]);
	if  length(c)>=23 & ~feof(fid1);  % a "valid" data line; otherwise, skip it
		if feof(fid1) | all( c(22:23)=='30');
			k2=0;
		else
			ste=c(22:23);
			%disp(['ste = ' ste]);
			if all(ste=='02') | all(ste=='29');
				disp(c(1:23));
				fprintf(fid2,'%s\n',c);
			end
		end
	else;
	end; % of if all(isstr(c)...
end


fclose(fid1);
fclose(fid2);

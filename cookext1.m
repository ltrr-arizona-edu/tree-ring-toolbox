function cookext1(i1)
% cookext1:  extract ascii subset of data for gridpoints from Ed Cook's pdsi
% CALL: cookext1(i1);
%
% Why?  Ed Cook's 154-gridpoint PDSI reconstruction.  Someone asked for the
% data for just a few gridpoints.  
%
%************** IN *************************
%
% i1 (1 x ?)i the desired gridpoints (e.g., [12 13 14]).
%
% User prompted for name of file (typically pdsirec1.dat) holding the pdsi
% reconstruction.  You should know the directory the file is in before 
% starting the function.  
%
%******************** OUT ********************
%
% No output arguments.
%
% User prompted for name of output file with the subset of data
% 



% **************** GET INPUT FILE

[file1,path1]=uigetfile('pdsirec*.dat','Infile with reconstructed PDSI');
pf1=[path1 file1];
fid1=fopen(pf1,'r');


% Note that gridpoint is right justified in cols 1921 of first header line
% and the reconstruction period is in cols 72-75 and 77-80 of next line
% Note also that 4 header lines for each gridpoint


% **************** GET OUTPUT FILE

[file2,path2]=uiputfile('*.dat','Outfile with subset reconstructed PDSI');
pf2=[path2 file2];
fid2=fopen(pf2,'w');





%************** CHECK INPUT

npnt = length(i1); % number of desired gridpoints

% Check that gridpoints in increasing order
diff1 = diff(i1);
if ~all(diff1>0);
	error('Specified gridpoints (i1) must be in ascending order');
end

% Check that gridpoints in range
if any(i1<1) | any(i1>155);
	error ('Gridpoints out of range');
end





%---------------   LOOP OVER GRIDPOINTS  

ntrack=0;
kstop=0;
while kstop==0;
	if ntrack==npnt;
		fclose all
		kstop=1;
		break
	end
	c=fgetl(fid1);
	if strcmp(c(18),'#'); % a header line 1 found
		ipoint = str2num(c(19:21)) % the gridpoint
		if ipoint == i1(ntrack+1);  % desired gridpoint found
         ntrack=ntrack+1;
         c2=fgetl(fid1);
			%Compute number of lines to copy
			yrgo = str2num(c2(72:75));
         yrsp = str2num(c2(77:80));
         dec1 = floor(yrgo/10);
         dec2= floor(yrsp/10);
         nln1 = dec2-dec1+1;
         nlines = 2 + nln1;
         
         fprintf(fid2,'%s\n',c);
         fprintf(fid2,'%s\n',c2);
			for n1 = 1:nlines;
				c=fgetl(fid1);
				fprintf(fid2,'%s\n',c);
			end
		end; % desired gridppoint not found
	end
			
end ; % of while
fclose all

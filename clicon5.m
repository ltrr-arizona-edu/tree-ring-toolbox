function clicon5
% clicon5: convert NCDC monthly "clivis" climate data to "decade line" format
%
% Meko 5-14-97
%
%
%**************  INPUT ****************************
%
%
% User is prompted for the following:
%
%  flnms*.txt -- list of files with the clivis data.  First line 
%			give source directory.  Second line gives target directory for the
%			output .mat files.  Remaining lines list the file suffixes.  
%
%			Example:
%
%				d:\greggin\
%				d:\greggout\
%				div07pds.dat
%				div07z.dat
%           div07tmp.dat
%           div07phi.dat
%
%				...etc
%
%
% The files div07pds.dat, etc., are ascii files in a weird format as
% downloaded from ncdc.   See notes
%
%
%************************ OUT *****************************
%
% * Set of .mat files, one per station, with year and monthly data in X
%   Data has these properties:
%
%		-13 cols (a year and jan-dec values) per line
%     -output .mat file with same prefix as input .dat file
%		-any missing gaps of years have been replaced with the year and 12 NaN
%		-missing monthly values are NaN
%
%
%********************* NOTES *********************************
%
% CLIVIS file format. 
%
% Line 1 -- blank
% Line 2 -- state  (e.g., Arizona)
% Line 3 -- blank
% Line 4 -- climatic division (e.g., Division 07)
% Line 5 -- blank
% Line 6 -- data type (e.g., Palmer Drought Severity Index)
% Line 7 -- blank
% Line 8 -- year (e.g., 1895)
% Lines 9-12 blank or unneeded title info
% Lines 13-24  -- paired month and value (e.g., 1  -4.59)
% Lines 25-26 -- blank
% Lines 27-43 -- sequence from year thru data
%
% repeat above thru all years
% 



%-------------- GET LIST OF FILES

[file1,path1]=uigetfile('flnms*.txt','Input list of filenames');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

% Get source and target paths from lines 1 and 2
c=fgetl(fid1);
if c(2)~=':';
	fclose all
	error('Line 1 is not valid path');
else
	path2=strtok(c); % path to input monthly data files
end

c=fgetl(fid1);
if c(2)~=':';
	fclose all
	error('Line 2 is not valid path');
else
	path3=strtok(c); % path to output .mat files
end


%-----------------  INTIAL READ TO COUNT NUMBER OF MONTHLY DATA FILES

nfiles=0;
k1=1;
while k1;
	c=fgetl(fid1);
	if ~feof(fid1)
		if any(c=='.'); % probably a valid file name
			nfiles =nfiles+1;
		else; % might be and eof on an empty line
			k1=0;
		end
	else;  % reached an eof
		k1=0;
	end
end
clc
disp(['Number of data filenames counted = ' int2str(nfiles)]);


% rewind and position at first filename
frewind(fid1);
c=fgetl(fid1);
c=fgetl(fid1);


%------- Initialize

% String matrix to hold file names
S1=blanks(12);
S2=repmat(S1,nfiles,1);




%************************* LOOP OVER INPUT DATA FILES

for n=1:nfiles
   
   % Get filename and build path filename
   cname=fgetl(fid1);
   pf2=[path2 cname];
   
   % Store filename prefix in a row of S2
   fpref=strtok(cname,'.');
   len1=length(fpref);
   S2(n,1:len1)=fpref;
   
   % Open the data file
   fid2=fopen(pf2,'r');
   if fid2==-1;
      fclose('all');
      error(['Cannot open ' pf2]);
   end
   
   
      
   % Report the state and the division to screen
   skipln1(fid2,1);
   cstate=fgetl(fid2);
   skipln1(fid2,1);
   cdiv = fgetl(fid2);
   skipln1(fid2,1);
   dtype=deblank(fgetl(fid2));
   skipln1(fid2,1); % should now be positioned for first year
   disp(['Working on ' cname ' -- ' cstate ' -- ' cdiv ]);
   
   
   
   % Initialize to hold 200 years of data -- more than enough
   a=NaN;
   X=repmat(a,200,13);
   
   
   
   
   % Open the output file
   pf3=[path3 fpref '.mat'];
   fid3=fopen(pf3,'w');
   
   % Loop over years
   nyears=0;
   k1=1;
   while k1;
      cyr = fgetl(fid2);
      if ~feof(fid2);
         year = str2num(cyr);
         skipln1(fid2,4);
         cdata=fscanf(fid2,'%g',24);
         X(nyears+1,:)=[year (cdata(2:2:24))'];
         nyears=nyears+1;
         skipln1(fid2,3);
      else
         k1=0;
      end
   end
   
   % Trim trailing all NaN years
   L1=(all(isnan(X')))';
   X(L1,:)=[];
   
   % Replace missing value code with NaN
   switch dtype
   case {'Palmer Drought Severity Index','Palmer Hydrological Drought Index'...,
      'Modified Palmer Drought Severity Index','Palmer Z Index'};
      mcode=999.99;
   case 'Temperature';
      mcode=-99.90;
   case 'Precipitation';
      mcode=-9.99;
   otherwise
      error('Invalid dtype');
   end
   
   L2=X==mcode;
   X(L2)=a;
   
      
   
   
   % Save output file
   eval(['save ' pf3 ' X']);
     
   fclose (fid2);
   fclose (fid3);
end
fclose(fid1);

   
   
   
   
   
   




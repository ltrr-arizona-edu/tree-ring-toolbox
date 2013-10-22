function rwlinp(pf1,path1,path4)
% rwlinp:  read .rwl-file ring widths and store in indexed-vector format in a .mat file
% rwlinp(pf1,path1,path4);
% Last revised 1-24-01 to handle rwl with series beginning before year -999 
%
% This is a way to get ring-width data into Matlab if the data is in the
% notorious ".rwl" format.  Index-vector form is simply each ring-width series
% stacked one after another in a single column vector.<P>
%
% Also stored is a matrix with years and row indices that allows retrieval of
% individual ring-width series from the I-V, and a string matrix of core ids
% in the same order as the ring-width series are stored. <P>
%
% Indexed-vector (I-V) storage is expected by several other functions in the
% tree-ring toolbox. An example is grplot.m, a function
% for plotting multiple ring-width series on a page.
%
%*** INPUT 
%
% pf1 (1 x ?)s path\filename of .rwl file
%   example:  'd:\jack\data\az033.rwl'
% path1 (1 x ?)s path to the .rwl file.  If this arg is passed, means
%   that pf1 is the filename only
%   Example: 'd:\jack\data\'  as path1,  and   'az033.rwl' as pf1
% path4 (1 x ?)s <optional, only if also have pf1 and path1>: path
%   for the output .mat files and .tmp files.  If no path4 as argument, 
%   default is to the same directory as the .rwl files are in
%
% Input arguments are optional.  There can be 0,1, 2, or 3 input arguments:
%  None: user prompted to clck on names of input and output files
%  One:  path\filename for the .rwl file
%  Two: first arg is the filename of the .rwl file, and the second is the path
%  Three: path for the output .mat file; this option is convenienent when user
%    wants .mat output files to go to different directories that that
%    of the source .rwl files
%
%*** OUTPUT  *****************************************
%
% No arguments. Depending on the number of input arguments, the output .mat
% file goes to a specified path\filename or the user is prompted to enter
% it.  The .mat file contains three variables:
%
%   X (? x 1)i  column vector of ringwidths stored one core after another
%	  in units of hundredths of mm
%   yrs (? x 3)i start year, end year, and row index of start year of each
%	  core's ring-width series in X
%   nms (? x 8)s  identification of each core
%
% A .tmp ascii file is also produced listing the core id, and first and last years'
% data for each core.  This file intended for checking that rwlinp.m
% indeed stores the ring widths properly
% 
%*** REFERENCES --none
%
%*** UW FUNCTIONS CALLED ******************
%
% intnan.m  -- checks for internal NaNs in a vector
% trailnan.m -- lops off trailing NaNs from a vector
%
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES ***************************************
%
% Handling of missing values in the .rwl files. Some ITRDB .rwl files
% have blanks for one or more internal years of a time series. In other
% words, the series is not continuous.  rwlinp.m checks for this, and
% warns the user, listing the series name and the violating year(s) on
% the screen.  The output X file will have NaNs in place of such 
% missing data.
%
% The core id.  rwlinp.m assumes the core id is in cols 1-8 of every line.
% A new series begins when any character of the the core id changes.
%
% Finding the first year of a ring-width series. rwlinp assumes the
% convention that the "year" on the first decade line of a ring-width
% series corresponds to the first year of data.  Like:
% 
% xxxxxx  1867   232   121    45
% xxxxxx  1870   etc    etc  etc  ...
%
% Note the start year above is 1867, and the three values on the first
% line are for 1867,68,69
%
% The end year of a ring-width series.  Assumed to be the last value
% preceding a data field with either all blanks or  a numeric 999
%
% End of file.  Assumed when any one of the following occurs:
%  1. eof encountered
%  2. ctrl-z starts a line
%  3. col 1 of a line is blank
% Some .rwl files have eof after the last data line. Some have eof on the
% last data line.  This is tested for in the initial count of data lines.
%
% Header lines in input .rwl.  ITRDB convention puts 3 header lines
% on the .rwl files.  rwlinp.m handles .rwl files with or without
% the header lines.  Header lines are assumed if the following 
% two conditions are satisfied:
%   1. col 8 of the first line in the file contains "1", col
%        8 of line 2 a "2", and col 8 of line 3 a "3"
%   2. lines 9:12 of the first line do not comprise a number
%
% A check is made for weird content of first data line.  I found that some
% non-ITRDB .rwl files have miscellaneous data in first line.  rwlinp
% checks that cols 9-12 of first data line 
% are consistent with being a year.  Col 9-12 must contain a single 
% right justified numeric value. Otherwise, error message. 
%
% Cols 1-8 generally contain a series code, like "CPE20A1 ", but series starting more than
% 1000 years before AD require 5 cols for the year (e.g., -1500).  Then the year is in cols 1-7.
%
% rwlinp.m is time consuming, but really only needs to be run once for a site.  Then
% the saved .mat file with data in I-V format can be loaded whenever you
% want to access ring widths for use in Matlab.
%
% See special calling program rwlinp1.m (not in toolbox) 
% to convert lots of .rwl files at one time
%
%------ assumed form of input .rwl file ---------------------
%	Col 1-7 or 1-8:  Core ID strings. 
%	Col 9-12 or 8-12 Year (decade)
%	Col 13-18, 19-24, ... Ring-width values, 10 per line
%
% The rwl file might optionally have a 3-line ITRDB header.  Function finds this
% out for checking that col 8 of lines 1,2,3 contain "1", "2", "3" respectively; 
% and that cols 9-12 of line 1 are inconsistent with being a year
%****************************************************************

% Initiate fit history
Fwhen = cell(8,4);
Fwhen{1,1}='rwlinp'; % function



%***********    GET INPUT .RWL PATH AND FILE NAME
if nargin==0; % interactive version; Prompt for data file name
   [fln,path1]=uigetfile('*.rwl','INPUT .rwl file');
   pf1=[path1 fln];
   path4=path1; % output to same directory as input .rwl files
else;  %automatic version;  path\file as input argument
   fslash=findstr(pf1,'\');
   if isempty(fslash); % no path prefixing filename; assume have 2 more args
      fln=pf1;
      if nargin~=3
 
         error(' Need paths to input .rwl and output .mat as args 2 and 3')
       end
   else; % need to extract input path from input argument pf1
      path1=pf1(1:max(fslash)); % overrides the path1 argument, if it exists 
      fln=pf1((max(fslash)+1):length(pf1));
      if nargin==1; % Output will go to same directory as input
         path4=path1;
      else; % Need 3 input args, and 3rd specifies dir for output
         if nargin ~=3,
      
            error('Need 3 input args for this setup -- third give out directory');
         end
      end
  end  
end

% Store name of input rwl file in fit history
Fwhen{1,3}=fln; % infile

% Open the data file for reading
fid=fopen(pf1,'r');

disp(['Input file: ' pf1]);

% Initialize
a=NaN;
blnks7=blanks(7);
blnks8=blanks(8); 
% Assume max possible 200 cores in a single rwl file
Nlines=a(ones(200,1),:);


% Find out if an ITRDB header  (3 lines) in file
c=fgetl(fid);
if length(c)<12;
	error('First line in .rwl file ends before col 12')
end
if c(8)=='1'  &   isempty(str2num(c(9:12)));
	ihead=1; % there is a set of 3 header lines
   c=fgetl(fid); % read 2nd header line
   if c(8)~='2';
      error('No 2 in col 8 of second header line');
   end
   c=fgetl(fid); % read 3rd headr line
   if c(8)~='3';
      error('No 3 in col 8 of third header line');
   end
   % now positioned for first data line
else
   yrcheck=str2num(c(9:12));
   if isempty(yrcheck) | length(yrcheck)~=1 | isspace(c(12))
      error('Cols 9-12 of first data line not a year')
   end
   
	ihead=0; % no header lines
	frewind(fid); % rewind to start of file
end

% Now positioned at first data line, either line 1 or lin 4 depending on ihead


%****************  READ PASS THRU DATA TO DETERMINE HOW MANY SERIES, AND 
% START AND END ROW OF EACH IN THE DATA FILE

disp('First pass: counting ring-width series and allocating space');

% start and end col for data values
jgo=(13:6:67)'; % start col for the 10 data values
jsp=jgo+5;
J1=[jgo jsp];  % start and end col

%idold=blnks8;
idold=blnks7;
k1=1;  % control on while loop
nsers=0; % series counter
linec=0; % line count
while k1;
	c1=fgetl(fid);
	% Check for first character being ctrl-z, or for end of file
   if feof(fid)==1 |  length(c1)<18 | isspace(c1(1)) ...
         | isempty(str2num(c1(9:12))) 
      if feof(fid)==1 & length(c1)>=12 & ~isempty(str2num(c1(9:12))); % Graybill type rwl
         % files have eof on end of last data line. Want to count that line.
         Nlines(nsers)=linec+1;
      else; % Regular .rwl files have eof after (below) last data line. Do
         % not want to count that line as data
         Nlines(nsers)=linec;
      end
      k1=0;
   else
      %idnew=c1(1:8);
      idnew=c1(1:7);
      if all(idnew==idold);  % same series
         linec=linec+1;
      else % first line of a new series
         if nsers==0;
            % Do not store count; this is first line of first series
            nsers=nsers+1;
            linec=1;
            idold=idnew;
         else ; % first line of series other than first series
				Nlines(nsers)=linec;  % store linecount for prev series
         	nsers=nsers+1;
            linec=1; % new series
            idold=idnew;
         end
         
      end
   end; % of if feof
end; % of while

% Trim off unneeded rows of Nlines
Nlines=Nlines(1:nsers);

% Using the total number of lines in the file, and considering 10 years max per
% line, compute an initializing length for the strungout vector
nmax = 10*sum(Nlines);

% Now know this:
%
% nsers is number of series
% Nlines(i) is the number of lines for each series
% Whether or no 3 header lines (ihead==1 vs ihead==0)

% Allocate
nms=blnks8(ones(nsers,1),:); % core ids
I=a(ones(nsers,1),:);
X=a(ones(nmax,1),:);   % strung out vector
yrs=a(ones(nsers,1),ones(3,1)); % start yr, end yr, row index of start year

%****************** REPOSITION TO START OF DATA IN FILE

frewind(fid);
if ihead==1;
	c=fgetl(fid);
	c=fgetl(fid);
	c=fgetl(fid);
end

%***** Store mask
cmask=ones(nsers,1);


%*************** LOOP OVER SERIES  *****************************


disp('Second pass: processing the data');

i=1; % target row in X
for n=1:nsers;
    if n==20;
        disp('20');
    end;
    
	%disp(['nsers = ' int2str(n)]);
	nrows=Nlines(n) ;  % number of decade lines for this core
	
	I(n)=i; % keep track of row in X of first value for this series
	% Loop over years
	for m = 1:nrows;
      c=fgetl(fid);
%       if n==41;
%           disp(['n = ' int2str(n)]);
%       end;
     
         
		% Handle first line
		if m==1;
            nmtemp = c(1:8);
            if nmtemp(8)=='-';
                nmtemp(8)=' ';
                yrslot1=8;
            else;
                yrslot1=9;
            end;
            nms(n,:)=nmtemp;
            yrgo=str2num(c(yrslot1:12));
            
            
            % Compute number years in initial decade
            if yrgo>=0;  % If 0 or A.D.
                ngo=10-rem(yrgo,10);
            else;
                ngo=-rem(yrgo,10);
                if ngo==0;  % A year like -310 as first year of decade should give 10 years in line
                    ngo=10;
                end;
                
            end
         
			% Read those years
			for j = 1:ngo;
				jgo=J1(j,1);
				jsp=J1(j,2);
				X(i)=str2num(c(jgo:jsp));
				i=i+1;
			end
		elseif m<nrows; % handle full decades
			for j=1:10,
				jgo=J1(j,1);
				jsp=J1(j,2);
				X(i)=str2num(c(jgo:jsp));
				i=i+1;
			end		
		else; % last row. either 10 values, or ends in 999 or  blanks
			decsp = str2num(c(yrslot1:12));   % decade 
			yrsp=decsp-1;  % tentative end year
			k2=1;
			j=1;
			while k2;
				jgo=J1(j,1);
				jsp=J1(j,2);
				if all(isspace(c(jgo:jsp))) | str2num(c(jgo:jsp))==999 | j>=10;
					k2=0;
				else
					X(i)=str2num(c(jgo:jsp));
					j=j+1;
					i=i+1;
					yrsp=yrsp+1;
				end
			end ; % of while k2
		end; % of if m== 
	end; % of for m= 
	yrs(n,:)=[yrgo yrsp I(n)];
end; % of for n over series

					
fclose(fid);


%*********************  CHECK CONSISTENCY OF YEAR RANGES AND ROW COUNTS

I=[I;i]; % running index of start row of series in target vector X
d1=diff(I); % should equal number of years in each series
d2=yrs(:,2)-yrs(:,1)+1; % likewise
d3=d2-d1;
if ~all(d3==0);
	error('Number of years of data inconsistent')
end
	
%************** Warning for leading or internal missing values in X

if isnan(X(1));
	disp('Leading NaN(s) in the strung out vector X');
	disp('Press any key to continue')
	disp(' ')
end

Lint=intnan(X);
if Lint==1; % trouble
	for n = 1:nsers;
		igo=yrs(n,3);
		isp=igo+yrs(n,2)-yrs(n,1);
		yr = (yrs(n,1):yrs(n,2))';
		x=X(igo:isp);
		f=find(isnan(x));
		if ~isempty(f); % here is trouble
			disp(['Core # ' int2str(n) '(' nms(n,:) ')']);
			disp('   Internal NaN in above series at years:');
			disp(yr(f));
			disp(' ');
			disp('Press any key to continue')
			pause
			disp(' ');
		end
	end
end
	

%**************** Get rid of any trailing NaNs from X ********

X=trailnan(X);

% Check that length of X consistent with row indices and years
% for last series

ilast = yrs(nsers,3)+yrs(nsers,2)-yrs(nsers,1);
if ilast ~= length(X);
	disp(['yrs vector says X should have row size ' int2str(ilast)]);
	disp(['  X has row size ' int2str(length(X))]);
	error('Row size of X inconsistent with yrs info (see above)');
end


%************************ MAKE OUTPUT FILE *********************

if nargin==0; % interactive version -- prompt for output file name
   txtfile=['.mat output file for ' fln]; 
   [file2,path2]=uiputfile('*.mat',txtfile);
   path4=path2;
else; % automatic verions -- assume .mat with same prefix as .rwl file,
   % and save to same directory as input file
   path2=path4;
   file2=strtok(fln,'.');
end

% Fit history
ctime=clock;
ctime=num2str(ctime(4:5));
dtime=date;
Fwhen{1,2}=[dtime ', ' ctime];
Fwhen{1,4}=file2;

% Save output
pf2=[path4 file2];
eval(['save ' pf2 ' X ' ' yrs ' ' nms cmask Fwhen']);


%*************************** MAKE SCREEN OUTPUT OF FIRST, LAST YEARS OFDATA

pf3=[path4 strtok(fln,'.') '.tmp'];
fid6=fopen(pf3,'w');

fprintf(fid6,'%s\n\n','First and Last Years of Ring-Width Data');
for n=1:nsers;
   which=nms(n,:);
   yrgo=yrs(n,1);
   yrsp=yrs(n,2);
   igo=yrs(n,3);
   isp=igo+yrsp-yrgo;
   xgo=X(igo);
   xsp=X(isp);
      
   str1=sprintf('%3.0f ',n);
   str2=sprintf('%s',which);
   str3=sprintf('%4.0f (%5.0f)   ',yrgo,xgo);
   str4=sprintf('%4.0f (%5.0f)',yrsp,xsp);
   strall=[str1 str2 str3 str4];
   
   fprintf(fid6,'%s\n',strall);   
end

fclose(fid6);


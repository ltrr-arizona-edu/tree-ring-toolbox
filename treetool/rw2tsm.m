function [X,yr,nms]=rw2tsm(path1,pf2,yrrange)
% rw2tsm: ring-width series to time series matrix
% [X,yr,nms]=rw2tsm(path1,pf2,yrrange);
% Last revised 7-19-99
%
% Reads ".rw files" of ring width from two or more cores and puts the ring-width
% series in a time series matrix (tsm).  The tsm by default covers the years
% from first year of earliest core through last year of most recent core.  The user 
% can optionally specify that the tsm is to cover some specific time period. 
%  
% The intended use is as a utility function to input the weirdly-formatted
% .rw ring-width data into Matlab for subsequent analysis.  The output arguments
% include series labels that can be used to identify series (e.g., in labeling plots).
%
% You must prepare an ascii file with a list of .rw filenames before running
% rw2tsm.  The tsm will have the ring-width series in the same order as the list
% of .rw files.
%
%*** IN ***********************
%
% path1 (1 x ?)s  path to .rw files  (e.g., c:\projs\ai4\treedata\)
% pf2 (1 x ?)s path\filename of file holding names of .rw files (one name pre line)
% yrrange(1 x 2)i  <optional> first and last year of desired time series matrix
%
%*** OUT ************************
%
% X (mX x nX)r ring widths for the nX series; each column a different series
% yr (mX x 1)i  year vector for X
% nms (nX x ?)s  string matrix of prefixes of .rw files. In same order as
%   columns of X
%
%*** REFERENCES-- none
%*** UW FUNCTIONS CALLED 
% rwread3 -- read individual .rw file
%
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES 
%
% If yrrange is specified, .rw data will be truncated to include only years in the 
% range.  An error is flagged if a .rw file has no data in the range.
%
% If yrrange is not specified, the output matrix X covers the inclusive period of the 
% set of .rw files.  In other words, from the earliest first year to latest last year.
% In allocating storage, I assume an initial coverage of AD 1 to AD 2050 (see hard code),
% but error flagged if a ring-width series falls outside the specified range
%
% NaN for missing data
% 
% The file pf2 can be written with notepad or another text editor. Be sure there are
% no trailing "blank" lines after the name of the last .rw file.  Be sure to include
% the suffix (".rw") in the file names.  
%
%************************

yrstart=1; yrend=2050;  % for initial allocation; needed if nargin==2

% Initialize storage of names
nms= 'blankline';

%************ Get file of .rw file names


% Read names of files to be copiedinto a column-cell of strings
s1 = textread(pf2,'%s','delimiter','\n','whitespace','');
%s1 = caseread(pf2);
[numrw,nchar]=size(s1);  % numrw is number of .rw files

% Put .rw file names into a cell array
s2=cellstr(s1); % cell array of .rw filenames

%******** Set number of rows (years) in output tsm
if nargin==3; % you have specified the year range for X
   nrow = yrrange(2)-yrrange(1)+1;
   yr = (yrrange(1):yrrange(2))';
else; % you want X to include all available data
   yrrange=[];
   yr= (yrstart:yrend)';
   nrow= length(yr);
end

X = repmat(NaN,nrow,numrw); % allocate


%  Loop over .rw files
for n = 1:numrw;
   file1 = s2{n};
   file1 = deblank(file1); % strip trailing blanks
   file1 =    fliplr(deblank(fliplr(file1))); % strip any leading blanks
   % Delete blank spaces
   L1 = isspace(file1);
   file1(L1)=[];
   % Store name
   nms=char(nms,strtok(file1,'.')); % store .rw file prefix 
      
   % put tsm of rw file into matrix X1, with year in col 1
   [X1,guy,day,fln]=rwread3(path1,file1);
   yrX1 = X1(:,1);
   yron=min(yrX1);
   yroff=max(yrX1);
        
   if nargin==3;% If yrrange was read in, put time series in proper col of X
      % First, last years of segments to pull
      yrgo = max([yron yrrange(1)]);
      yrsp = min([yroff yrrange(2)]);
      
      % Pointers
      L1 = yrX1>=yrgo & yrX1<=yrsp;  % to X1
      L2 = yr>=yrgo & yr<=yrsp; % to X
      if ~any(L1);
         error(['Series ' int2str(n) ' not in desired year range']);
      end
      
      % Pull data
      x=X1(L1,2);
      
      % Store data
      X(L2,n)=x;
      
   else; % want X to hold all available ring widths
      if yron < yrstart;
         error(['Too early a start year for series ' file1]);
      end
      if yroff > yrend;
         error(['Too late an ending year for series ' file1]);
      end
      L1 = yr>=yron & yr<=yroff;
      X(L1,n)=X1(:,2);
   end
end; % of loop over .rw files

% Lop off all-NaN rows
L = (all(isnan(X')))';
if any(L);
   X(L,:)=[];
   yr(L)=[];
end

% Delete leading row of nms
nms(1,:)=[];
      
      

function [x,yrx,n]=getts(flnm,datatype,tsname)
% getts: get a time series from meko storage file
% CALL: [x,yrx,n]=getts(flnm,datatype,tsname);
% 
% Meko 6-28-98
%
%********* IN
%
% flnm (1 x ?)s name (w/0 suffix) of mat file storing chronology
%   (e.g., c:\projs\ah4\treedata\rw\med)
% datatype (1 x ?)s data type
%    X -- ring width
%    IX -- core index, standard
%    EX -- core index, residual
%    IT -- tree index, standard
%    ET -- tree index, residual
%    ZI -- site index, standard
%    ZE -- site index, residual
% tsname (1 x ?)s name of series;  content will vary with datatype
%   ZI or ZE: empty [] 
%   X, IX, or EX: core name as used in nms variable in storage file
%   IT or ET:  tree name as used in Tnms variable in storage file
%
%********** OUT 
%
% x (1 x ?)r  core index
% yrx (1 x ?)i corresponding year vector
% n (1 x ?)i number of "samples" in each year of x;  This meaningful only for
%   datatypes IT,ET, ZI or ZE.  For tree indices, indicates number of cores
%   making up the tree index. For site chron, indicates number of trees making
%   up the site index
%
%********* NOTES
%
% Check beforehand to make sure you know the series name spelling as
% as used in the .mat file



%-----  LOAD THE .MAT FILE

% Strip .mat off filename if it is on there
if ~isempty(findstr(flnm,'.'));
   flnm=strtok(flnm,'.');
end

flnm = [flnm '.mat'];
if exist(flnm)~=2;  
   error('Requested .mat storage file flnm does not exist');
end
eval(['load ' flnm]);

%-------- GET NEEDED INFORMATION

% Ringwidth
if strcmp(datatype,'X');
   row1 = strmatch(tsname,nms);
   if isempty(row1);
      error(['no ' tsname ' in nms']);
   end
   nfound = length(row1);
   if nfound>1;
      error([int2str(nfound) ' series match name ' tsname]);
   end
   % Pull row index information for series
   years = S(row1,10:11); % start and end year of series
   igo = S(row1,12); % row index of first year in storage matrix
   isp = igo + diff(years); % computed row index of last year
   x = X (igo:isp);  % desired data
   yrx = (years(1):years(2))'; % year vector
   n=[];
elseif any(strcmp(datatype,{'IX','EX'})); % want standard or residual core index
   yrs=yrs; % years and row index info, same for both data types
   if strcmp(datatype,'IX');
      X = 'IX'; % data matrix, standard core index
   else
      X = 'EX'; % data matrix, residual core index
   end
   row1 = strmatch(tsname,nms);
   if isempty(row1);
      error(['no ' tsname ' in nms']);
   end
   nfound = length(row1);
   if nfound>1;
      error([int2str(nfound) ' series match name ' tsname]);
   end
   % Pull row index information for series
   years = yrs(row1,1:2); % start and end year of series
   igo = yrs(row1,3); % row index of first year in storage matrix
   isp = igo + diff(years); % computed row index of last year
   eval(['x = ' X '(igo:isp);'])  % desired data
   yrx = (years(1):years(2))'; % year vector
   n=[]
   
elseif any(strcmp(datatype,{'IT','ET'})); % want standard or residual tree index
   if strcmp(datatype,'IT');
      X = 'IX'; % data matrix, standard tree index
      yrs = 'ITyrs' % years and row info
      XN= 'ITn'; % number of cores in tree
   else
      X = 'EX'; % data matrix, residual core index
      yrs = 'ETyrs'; % years and row info
      XN = 'ETn'; % number of cores in tree
   end
   row1 = strmatch(tsname,Tnms);
   if isempty(row1);
      error(['no ' tsname ' in Tnms']);
   end
   nfound = length(row1);
   if nfound>1;
      error([int2str(nfound) ' series match name ' tsname]);
   end
   % Pull row index information for series
   eval(['years = ' yrs '(row1,1:2);']); % start and end year of series
   eval(['igo = ' yrs '(row1,3);']); % row index of first year in storage matrix
   isp = igo + diff(years); % computed row index of last year
   eval(['x = ' X '(igo:isp);'])  % desired data
   eval(['n = ' XN '(igo:isp);'])  % desired info on number of cores in each yr
   yrx = (years(1):years(2))'; % year vector
   
   
elseif any(strcmp(datatype,{'ZI','ZE'})); % want standard or residual site index
   if ~isempty(tsname);
      error('tsname should be empty in calling for ZI or ZE');
   end
   if strcmp(datatype,'ZI');
      X = 'ZI'; % data matrix, standard site index
      yrs = 'yrZI' % year vector
   else
      X = 'ZE'; % data matrix, residual site index
      yrs = 'yrZE'; % years and row info
   end
   
   eval(['yrx = ' yrs ';']); % year vector
   eval(['x = ' X '(:,1);']); % data
   eval(['n = ' X '(:,6);']); % number of trees in site chron each year
      
else
   error(['Datatype ' datatype ' invalid']);
end

   



 
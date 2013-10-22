function rwc2rw(XT,XEL,path1,ftype,kopt)
% rwc2rw: ring-width cells to .rw, .eww or .lww files
% rwc2rw(XT,XEL,path1,ftype);
% Last revised 2011-04-04
%
% Converts cell-format ring-width data to .rw, .eww or .lww files.  Called from edit mode of ring-width measuring 
% function rwmeas.m 
%
%
%*** INPUT
%
% XT, XEL structure variables of cells containing ring-width measurements and associate information (see rwmeas.m)
% path1 (1 x ?)s   path to write output files to
% ftype (1 x ?)s   type of file to create
%   == rw  .rw files of total measured ring width
%   == ew/lw  .eww and .lww files of measured earlywood and latewood 
%   == rwc  .rw files of computed total ring width (EWW+LWW)
%   == rwm  .rw files of merged total-ringwidth data from XT and XEL (see notes)
% kopt (1 x 1)i optional, tells whether want rw data as hundredths or thousandths of mm
%   ==1 hundredths  (also if no kopt -- nargin==4)
%   ==2 thousands
%
%*** OUTPUT
%
% No arguments. 
% .rw, .eww or .lww files are written
%
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED 
%
% trailnan -- removes trailing NaNs
% uniqcell -- finds unique series names
% writerw -- writes a .rw, .eww or .lww file
%
%
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% rwv option:  total ring width data exclusively taken from measured total ring width
% rwc option: total ring width data exclusively taken from computed (EW+LW) data (file XEL)
% rw option. total ring width data is taken from the partial-width data (EW+LW), or if that is not available, from the
%   measured total-width data (XT).  Thus the priority is for the partial-width measurements.  If data exist for both
%   total and partial width for a given year, priority is given to the partial-width data.
%
% Rev2006-08-09:  to handle thousandths of mm output.  Now multiplies measurements by fscale rather than by
% 100.
% Rev2011-04-04:  to handle path and file naming problems in linux
%

%--- UNITS OF OUTPUT RWL

if nargin<5;
    kopt=1;
else;
    L=kopt==[1 2];
    if ~any(L);
        error('kopt must be 1 or 2');
    end;
end
switch kopt;
    case 1;
        fscale=100;
    case 2;
        fscale=1000;
end



% check path1

if ~isstr(path1);
    error('path1 not a string');
end;
if size(path1,1)~=1;
    error('path1 has row dimension >1');
end;

% Rev2011-04-04:  to handle path and file naming problems in linux
% % if ~isletter(path1(1));
% %     error('path1 must begin with a letter');
% % end;
% % if ~strcmp(path1(end),'\');
% %     error('path1 must end with a backslash');
% % end;
if ~(isletter(path1(1))  ||  strcmp(path1(1),'/'));
    error('path1 must begin with a letter or forward slash');
end;
if ~(strcmp(path1(end),'\')  ||  strcmp(path1(end),'/'));
    error('path1 must end with a backslash or forward slash');
end;





% Check XT, XEL

if ~isstruct(XT);
    error('XT must be a structure variable');
end
if ~all([isfield(XT,'data')  isfield(XT,'when')  isfield(XT,'who')  isfield(XT,'id')]);
    error('XT does not have the required fields');
end;

if ~isstruct(XEL);
    error('XEL must be a structure variable');
end
if ~all([isfield(XEL,'data')  isfield(XEL,'when')  isfield(XEL,'who')  isfield(XEL,'id')]);
    error('XEL does not have the required fields');
end;

%-- SET SUFFIX OF FILES

switch ftype;
case {'rw','rwc','rwm'}; %  total ringwidth
    suffx = 'rw';
case {'ew/lw'}; 
    suffx1='eww';
    suffx2='lww';
otherwise;
    error('ftype invalid -- must be rw, ew/lw or rwc');
end; % switch ftype


% -- FIELDS TO CELL VARIABLES
id1 = XT.id;
d1 = XT.data;
w1=XT.when;
p1= XT.who;
id2 = XEL.id;
d2 = XEL.data;
w2=XEL.when;
p2= XEL.who;


% --- CHECK THAT AT LEAST ONE SERIES HAS DATA FOR TOTAL OR EW/LW

switch ftype;
case {'ew/lw','rwc'};
    nser = length(id2);
    if nser<=1;
        if isempty(id2{1});
            error('No series measured yet for EWW/LWW (XEL empty)');
        end;
    end;
case {'rwm'};
    nser = length(id1);
    if nser<=1;
        if isempty(id1{1});
            error('No series measured files yet for total width (XT empty)');
        end;
    end;
case 'rw';
    nser = max([length(id1) length(id2)]);
    if nser<=1;
        if isempty(id1{1}) & isempty(id2{1});
            error('No measured files yet in XT or XEL');
        end;
    end;
    if isempty(id1{1}); % no total width data
        ftype='rwc'; % must use computed total width 
    elseif isempty(id2{1});
        ftype = 'rwm'; % must use measured only total
    end;
    
otherwise;
    error([ftype ' invalid ftype']);
end;



% -- CALL TO WRITE FILES

switch ftype;
case 'rwm';
    for n = 1:nser;
        X=d1{n};
        yrx=X(:,1);
        x = round(fscale*X(:,2));
        
        % Strip trailing and leading NaN
        [x,yrx]=subfcn01(x,yrx);
                
        yrfirst = yrx(1); % first year of valid data
        person = p1{n}; % measurer
        dwhen = w1{n}; % date measured, or last edited
        fn = lower(id1{n}); % filename set to id of series
        writerw(x,person,dwhen,yrfirst,fn,suffx,path1);
    end;
case 'rwc';
    for n = 1:nser;
        X=d2{n};
        yrx=X(:,1);
        x = round(fscale*X(:,4));
        
        % Strip trailing and leading NaN
        [x,yrx]=subfcn01(x,yrx);
                
        yrfirst = yrx(1); % first year of valid data
        person = p2{n}; % measurer
        dwhen = w2{n}; % date measured, or last edited
        fn = lower(id2{n}); % filename set to id of series
        writerw(x,person,dwhen,yrfirst,fn,suffx,path1);
    end;
case 'rw';
    
    idset = uniqcell(id1,id2); % get cell array of series ids, no duplicates, from XT and XEL
    nset = size(idset,2); % number of ids
    
    for n = 1:nset;
        
        %-- Get time series of total width, measured or computed
        idthis = idset{n};
        j1 = strcmpi(idthis,id1);
        if ~all(j1==0); % measured total width exists for this series
            X1=d1{j1}; % total width data
            yrx1 = X1(:,1);
            x1=round(fscale*X1(:,2));
            % Strip trailing and leading NaN
            [x1,yrx1]=subfcn01(x1,yrx1);
            % Set other series info
            person=p1{j1};
            dwhen=w1{j1};
            dwhen=w1{j1};
            fn = lower(id1{j1}); 
            
        else;
            x1=[];
        end;
        j2 = strcmpi(idthis,id2);
        if ~all(j2==0);
            X2=d2{j2}; % total width computed from partial
            yrx2=X2(:,1);
            x2=round(fscale*X2(:,4));
            % Strip trailing and leading NaN
            [x2,yrx2]=subfcn01(x2,yrx2);
            % Set other series info
            person=p2{j2};
            dwhen=w2{j2};
            dwhen=w2{j2};
            fn = lower(id2{j2}); 
        else;
            x2=[];
        end;
        
        
        if isempty(x1) & isempty(x2);
            error('Neither partial nor total data for an identified series');
        elseif isempty(x1); % if no total-width measurements, use the partial
            x = x2; 
            yrx=yrx2;
        elseif isempty(x2); % if no partial-width measurements, use the total
            x=x1;
            yrx=yrx1;
        else; % apparently have both types of data
            % Compute start year for merged series
            yrgo = min([yrx1(1) yrx2(1)]);
            yrsp = max([yrx1(end) yrx2(end)]);
            
            % Initialize x and yrx
            yrx = (yrgo:yrsp)';
            nyr = length(yrx);
            x = repmat(NaN,nyr,1);
            
            % Substiture measured total
            i1= yrx1-yrgo+1;
            x(i1)=x1;
            
            % Substitute computed total (EW+LW) into x, overwriting any measured total width for the same years
            i2 = yrx2-yrgo+1; % row indices
            x(i2)=x2;
        end; % if isempty(x1) & isempty(x2);
        
        yrfirst = yrx(1); % first year of valid data
        
        writerw(x,person,dwhen,yrfirst,fn,suffx,path1);
    end;
    
    
case 'ew/lw';
    
    % eww
    for n = 1:nser;
        X=d2{n};
        yrx=X(:,1);
        x = round(fscale*X(:,2)); % to hundredths of mm
        
        % Strip trailing and leading NaN
        [x,yrx]=subfcn01(x,yrx);
        
        yrfirst = yrx(1); % first year of valid data
        person = p2{n}; % measurer
        dwhen = w2{n}; % date measured, or last edited
        fn = lower(id2{n}); % filename set to id of series
        writerw(x,person,dwhen,yrfirst,fn,suffx1,path1);
    end;
    % lww
    for n = 1:nser;
        X=d2{n};
        yrx=X(:,1);
        x = fscale*round(X(:,3));
                
        % Strip trailing and leading NaN
        [x,yrx]=subfcn01(x,yrx);
        
        yrfirst = yrx(1); % first year of valid data
        person = p2{n}; % measurer
        dwhen = w2{n}; % date measured, or last edited
        fn = lower(id2{n}); % filename set to id of series
        writerw(x,person,dwhen,yrfirst,fn,suffx2,path1);
    end;
otherwise;
end;


disp(' ');
switch ftype;
case {'rw','rwc','rwm'}
    disp([int2str(n) ' .rw files successfully written to directory ' path1]);
case 'ew/lw';
    disp([int2str(n) ' .eww and ' int2str(n) ' .lww files successfully written to directory ' path1]);
end;
    

%---- SUBFUNCTION 01    
function [x,yrx]=subfcn01(x,yrx);
% Strip trailing and leading NaN
x=trailnan(x);
mx=length(x);
yrx=yrx(1:mx);
x=flipud(x);
yrx=flipud(yrx);
x=trailnan(x);
mx=length(x);
yrx=yrx(1:mx);
x=flipud(x);
yrx=flipud(yrx);




















    



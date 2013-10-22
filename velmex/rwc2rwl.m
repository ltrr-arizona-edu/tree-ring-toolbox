function rwc2rwl(XT,XEL,path1,ftype,rwlset,kopt);
% rwc2rwl:  cell-format ringwidth data to .rwl files
% rwc2rwl(XT,XE,pathout,ftype,rwlset,kopt);
% Last revised 2011-3-14
%
% Make .rwl files (decade format) of total ring width, earlywood width, or latewood width from cell-format 
% ring-width storage files (see rwmeas.m)
%
%*** INPUT
%
% XT, XEL structure variables of cells containing ring-width measurements and associate information (see rwmeas.m)
% path1 (1 x ?)s   path to write output files to
% ftype (1 x ?)s   type of .rwl file to create
%   == rw  total ringwidth, using all available data
%   == ew  earlywood width width 
%   == lw  latewood width 
%   == rwc  total ring width, computed only EWW+LWW)
%   == rwm  total ring width, measured directly only (not using any EW+LW)
% rwlset -- structure with information on tailored sets of series to be included in rwl files
% kopt (1 x 1)i  <optional>  units of desired output rwl measurement (see notes)
%   1==hundredths of mm
%   2==thousandths of mm
%   no arg --> just like kopt==1
%
%*** OUTPUT
%
% No arguments. 
% .rwl files are written
%
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED 
% tsv2rwl -- time series vector to block of rwl data
%
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% rw option:  total ring width 
%   ring width data of first priority is measured total ring width (col 2 in XT);  second priority is
%   computed total ring width (col 4 of XEL). Target vector first filled with the computed total ring width. Then
%   data overwritten with measured total ring width if measured total ring-width available for those years
%
% rwc option: total ring width data exclusively taken from computed (EW+LW) data (file XEL)
% rwm option. data exclusively taken from measured total ring width 
%
% kopt:  Revised 2006-8-9 for Scott St. George to include the thousandths of mm digit% of ring width, optionally. The older
% rwl files in the ITRDB are integers in hundredths of mm (e.g., 147 is 1.47 mm).  Machines now measure to
% thousandths of mm, and some newer ITRDB files have ringwidth to thousandths (e.g., 1473 is 1.473 mm). It is
% up to you to decide how you want your rwl files to represent data.  My Velmex stage measures to thousandths.
% If I set kopt=1 or do not include kopt as an argument in call to rwc2rwl, my internal measurement of 1.473 mm
% will print as 147.  If kopt=2, the measurement will print in the rwl as 1473.  
%
% Rev2011-3-15.    To accept a "/" as start of path. Otherwise trips up in Linux.  Previously demanded that path began as letter. 


% FIND OUT WHETHER HUNDREDTHS OR THOUSANDS OF MM AS OUTPUT

if nargin<6;
    kopt=1;
else
    L=kopt==[1 2];
    if ~any(L);
        error('kopt, if present as input arg, must be 1 or 2');
    end
end;
switch kopt;
    case 1;
        fscale=100;
    case 2; 
        fscale=1000;
    otherwise
end




% CHECK PATH

if ~isstr(path1);
    error('path1 not a string');
end;
if size(path1,1)~=1;
    error('path1 has row dimension >1');
end;

% Rev2011-3-15
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


%-- PROMPT FOR FIRST 3 CHARS OF OUTPUT FILENAME 

kwh2=1;
while kwh2==1;
    prompt={'Enter characters (first must be a letter):'};
    def={'XXX'};
    dlgTitle='First 3 characters of desired .rwl output filename';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    ch3 = answer{1};
    if ~all([length(ch3)==3    isletter(ch3(1))  ~any(isspace(ch3))   ]);
        uiwait(msgbox('Enter 3 characters again','Invalid ch3','modal'));
    else;
        kwh2=0;
    end;
end;




%-- SET SUFFIX OF FILES

suffx='rwl';


% -- FIELDS TO CELL VARIABLES
% id is series id
% w is date measurements completed
% p is person measuring
% d is the data (the measurements
% 1 refers to XT data (measured total ringwidth)
% 2 refers to the XEL data (partial width and computed total)

id1 = XT.id;
d1 = XT.data;
w1=XT.when;
p1= XT.who;
id2 = XEL.id;
d2 = XEL.data;
w2=XEL.when;
p2= XEL.who;



%----   PROMPT FOR OPTIONAL USE OF A TAILORED RWLSET; IF USING ONE, CULL THE SUBSET DATA

krwlset=menu('Choose','Use subset of series according to a previously tailored rwlset',...
    'Use all series');
if krwlset==1; % Use an rwlset
    [id1,d1,id2,d2]=cullids(id1,d1,id2,d2,rwlset,ftype);
else;
end;



% --- CHECK THAT AT LEAST ONE SERIES HAS DATA FOR TOTAL OR EW/LW

switch ftype;
case {'ew','lw','rwc'};
    nser = length(id2);
    if nser<=1;
        if isempty(id2{1});
            error('No measured files yet in XEL');
        end;
    end;
case {'rwm'};
    nser = length(id1);
    if nser<=1;
        if isempty(id1{1});
            error('No measured files yet in XT');
        end;
    end;
case 'rw';
    nser = max([length(id1) length(id2)]);
    if nser<=1;
        if isempty(id1{1}) & isempty(id2{1});
            error('No measured files yet in XT or XEL');
        end;
    end;
    % Number of series to treat is number of nonduplicate ids from id1 and id2
    idset = uniqcell(id1,id2); % get cell array of series ids, no duplicates, from XT and XEL
    nset = size(idset,2); % number of ids
    nser=nset;
    
otherwise;
    error([ftype ' invalid ftype']);
end;



% WRITE FILES 


% str1 == default ending part of file stem (e.g., padwe1.rwl has we1 indicating width, earlywood, version1

%-- STORE YEAR VECTOR AND DATA, TREATING DATA DIFFERENTLY DEPENDING ON FTYPE
switch ftype;
case 'rwm';
    str1='wt1';
    Xkey = d1;
    keycol=2;
        
case 'ew';
    str1='we1';
    Xkey = d2;
    keycol=2;
        
case 'lw';
    str1='wl1';
    Xkey = d2;
    keycol=3;
    
case 'rwc';
    str1='wt1';
    Xkey = d2;
    keycol=4;
      
case 'rw'; % Merge data from XT and XEL
    str1='wt1';
    Xkey = [];
    keycol=[];
    
end; % switch ftype 



%--  BUILD STRING MATRIX OF RWL-FORMATTED DATA

if strcmp(ftype,'ew') | strcmp(ftype,'lw') | strcmp(ftype,'rwc');
    for n = 1:nser; 
        nm = upper(id2{n}); % id of series
        
        X=Xkey{n};
        x = round(fscale*X(:,keycol));
        yrx=X(:,1);
        
        % Strip trailing and leading NaN
        [x,yrx]=subfcn01(x,yrx);
        
        % call tsv2rwl to get block of char data
        sthis=tsv2rwl(nm,x,yrx);
        
        % Paste new block onto existng 
        if n==1;
            S=sthis;
        else;
            S=char(S, sthis);
        end;
        
    end; % for n = 1:nser;
elseif strcmp(ftype,'rwm'); % total width, measured only
    for n = 1:nser; 
        nm = upper(id1{n}); % id of series
        
        X=Xkey{n};
        x = round(fscale*X(:,keycol));
        yrx=X(:,1);
        
        % Strip trailing and leading NaN
        [x,yrx]=subfcn01(x,yrx);
        
        % call tsv2rwl to get block of char data
        sthis=tsv2rwl(nm,x,yrx);
        
        % Paste new block onto existng 
        if n==1;
            S=sthis;
        else;
            S=char(S, sthis);
        end;
        
    end; % for n = 1:nser;
    
elseif strcmp(ftype,'rw'); % ftype == 'rw';
    
    for n=1:nser; % loop over series
        
        %-- Get time series of total width, measured or computed
        idthis= idset{n};
        j1 = strcmpi(idthis,id1);
        if ~all(j1==0); % measured total width exists for this series
            X1=d1{j1}; % total width data
            yrx1 = X1(:,1);
            x1=round(fscale*X1(:,2));
            % Strip trailing and leading NaN
            [x1,yrx1]=subfcn01(x1,yrx1);
            
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
        else;
            x2=[];
        end;
        
        if isempty(x1) & isempty(x2);
            error('Neither partial nor total data for an identified series');
        elseif isempty(x1); % if no total-width measurements, use the computed from partial
            x = x2; 
            yrx=yrx2;
        elseif isempty(x2); % if no partial-width measurements, use the measured-only total
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
            
            % Substitute measured total
            i1= yrx1-yrgo+1;
            x(i1)=x1;
            
            % Substitute computed total (EW+LW) into x, overwriting any measured total width for the same years
            i2 = yrx2-yrgo+1; % row indices
            x(i2)=x2;
        end; % if isempty(x1) & isempty(x2);
        
        % call tsv2rwl to get block of char data
        nm=idthis;
        if strcmp(nm,'PRC08A');
            disp('here');
        end;
        
        sthis=tsv2rwl(nm,x,yrx);
        
        % Paste new block onto existng 
        if n==1;
            S=sthis;
        else;
            S=char(S, sthis);
        end;
        
        
    end; % loop over n=1:nser
   
end; % if ~strcmp(ftype,'rw');


% SET OUTPUT FILENAME
deffn = [ch3 str1 '.' suffx]; % default filename
pfdef =[path1 deffn];
kmen=menu('Choose one',...
    ['Proceed writing ' pfdef],...
    'Select a different output .rwl filename');
if kmen==1;
    fnout=deffn;
else
    pause(2)
    [file2,path2]=uiputfile('*.rwl',['Output file instead of ' pfdef]);
    pfdef=[path2 file2];
end;


%-- WRITE OUTPUT FILE
fid1=fopen(pfdef,'w');
[mS,nS]=size(S);
fmt = '%s\r\n';
for n = 1:mS;
    s= deblank(S(n,:));
    fprintf(fid1,fmt,s);
end;
fclose(fid1);

%-- MESSAGE TO USER
uiwait(msgbox([pfdef ' sucessfully written'],'Message','modal'));


%---- SUBFUNCTION 01    ---STRIP TRAILING AND LEADING NAN FROM TIME SERIES
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



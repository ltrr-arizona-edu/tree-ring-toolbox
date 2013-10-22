function rwmeas
% rwmeas:  ring-width measurement
% rwmeas;
% Last revised 2-11-02
%
% Measure ring widths with Velmex system.  Options for earlywood/latewood width or total width
%
%*** INPUT
%
% No args.  User prompted for files and settings
%
%
%*** OUTPUT 
%
% No args. User prompted for filenames.  File created or modified is .mat storage of 3-d matrix
% with measured series.  
%
%
%*** NOTES
%
% Can append data to .mat file, or create new file
% Use other functions ?? for converting data to .rw or .rwl
% Data to be stored:
%   XT, XEL  -- structures 
%      Fields: each field is a structure variable 
%       .units -- string measurement units (e.g., 'mm x 100')
%       .summary -- char array listing series ids and periods covered
%       .data -- cell with tsm for each measured series, year as col 1. (2 col for total, 4 col for EL)
%           For example, XT.data{1} has a tsm of total width for series
%       .id --- cell with id of each series
%       .span -- first, last year of series 
%       .who -- measurer -- also a cell variable
%       .when -- date measurement completed --- also in cells
%       .remeasure -- years of data flagged to be re-measured (e.g., after editing in edit mode to split a ring
%  rwlset  --structure with specification for series/periods for rwl files; for jth rwl set:
%       .names{j} (1 x ?)s short name name (e.g., 12 chars) 
%       .describe{j} (1 x ?)s  text description of the rwl set
%       .when{j}   date the rwlset created or modified
%       .idnames{j} cell with ids for each series in the set
%       .trimall{j} (1 x 2)i  start and end year for the rwl set (used to limit the period for all series)
%       .trimeach{j} cell of (1x2)i with specified start and end year for each series in the rwl set
%  sinfo -- string matrix summarizing years of measurments for each series
%  vlist (?x?)s  definitions of stored variables
%
% For example, XEL.data{3} is a 4-column tsm of year , EWW, LWW, EWW+LWW for series 3. 
% XEL.span{3} is a  3x2 matrix with first and last years (span) of EWW, LWW, EWW+LWW 
% XEL.who{3} holds a string variable of who measured --- e.g.,   'DMM'
% XEL.when{3} holds date and time last measurements on the series made -- e.g., 08-05-01
% XEL.data{} has 4 columns: year, EWW, LWW and total of EWW+LWW.  This last column of total ring width contrasts with
%   the contents of XT, which is directly measured total ring width.
%
% Series are sorted and stored in alphabetical order by id

clear;
close all;
clc;

%--- PROGRAM SETTINGS (CAN HARD CODE THESE AS NEEDED)

prompt={'Enter your initials','Enter the units for measurments as stored in files',...
        'Enter maximum allowable series length (yr)',...
        'Enter number of leading NaN to allow extension of measurements before first existing measured year'};
def={'skn','mm','2000','200'};
dlgTitle='Program Settings';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
person=answer{1};
units=answer{2};
maxlen=str2num(answer{3});
nyrlead=str2num(answer{4});


%---- MEASURE OR EDIT?

kmen=menu('What doe you want to do?','Measure','Edit');
if kmen==1;
    runmode='Measure';
else;
    runmode='Edit';
end;
clear kmen;

switch runmode;
case 'Edit';  % EDIT MODE
    
    % Menu for editing desired
    kedit1=menu('Choose',...
        'Edit measurements for single series',...
        'Bulk edit multiple series',...
        'Import .rw, .eww or .lww files         ',...
        'Create or revise an rwl specification set',...
        'Write output data files'); 
    
    kfile='Existing'; % by default, work with an existing .mat cell-format ring width file
    % Loading the existing .mat file or making new one
    if kedit1==3; % if importing dat
        kmen6= menu('Choose one import mode',...
            'Import into an existing .mat storage file',...
            'Create a new .mat storage file for importing data to');
        if kmen6==1; % import into an existing
            % kfile remains 'Existing';
        else;
            kfile='New';
            
        end;
        
    end;
    filemode=kfile;
    
    % LOAD OR CREATE THE .MAT STORAGE FILE
    
    if ~strcmp(kfile,'New');
        [file1,path1]=uigetfile('*.mat','Input file storing ring-width data');
        pf1=[path1,file1];
        eval(['s=load(''' pf1 ''');']);
        if ~all([isfield(s,'vlist')  isfield(s,'XEL')  isfield(s,'XT')  isfield(s,'sinfo') isfield(s,'rwlset')]);
            error([pf1 ' does not contain XEL, XT , sinfo, rwlset and vlist']);
        end
        clear s;
        % Load the ring-width data
        eval(['load ' pf1 '  XT XEL sinfo rwlset vlist;']); 
    else;
        [XT,XEL,sinfo,rwlset,vlist,pf1]=subfcn08(units);
    end;
    
   
    
    % IF THROUGH THE EDIT CHOICES
        
    if kedit1==1; % Edit ring-width measurements
        filemode='Existing';  % Cannot edit a new file
        
        %--- TYPE OF DATA TO EDIT
        kmen1=menu('EDIT WHAT?',...
            'TOTAL WIDTH',...
            'PARTIAL WIDTH');
        if kmen1==1;
            editmode='Total';
            XT=rwedit01(XT);
            vbl='XT';
        else;
            editmode='Partial';
            XEL=rwedit01(XEL);
            vbl='XEL';
        end;
        
        clear kmen1;
    elseif kedit1==2; % bulk edit
        filemode='Existing'; % cannot bulk edit that which does not exist
        %--- TYPE OF DATA TO EDIT
        kmen1=menu('EDIT WHAT?',...
            'TOTAL WIDTH',...
            'PARTIAL WIDTH');
        if kmen1==1;
            editmode='Total';
            XT=editbulk(XT);
            vbl='XT';
        else;
            editmode='Partial';
            XEL=editbulk(XEL);
            vbl='XEL';
        end
           
    elseif kedit1==3;  % import data from a .rw, .eww or .lww file measured outside the system
        
        %--- TYPE OF DATA TO EDIT
        kmen1=menu('EDIT WHAT?',...
            'TOTAL WIDTH',...
            'PARTIAL WIDTH');
        if kmen1==1;
            editmode='Total';
            XT=rwimport(XT,editmode);
            vbl='XT';
        else;
            editmode='Partial';
            XEL=rwimport(XEL,editmode);
            vbl='XEL';
        end;
        
    elseif kedit1==4; % make or revise rwlset
        rwlset=rwlspecs(rwlset,XT.id,XEL.id);
        
        
    elseif kedit1==5;;  % create .rw, etc files 
        pathdef = path1;
        prompt={'Enter the path for output files:'};
        def={pathdef};
        dlgTitle='Path for output';
        lineNo=1;
        answer=inputdlg(prompt,dlgTitle,lineNo,def);
        pathout = answer{1};
        
        kwh1=1;
        while kwh1==1;
            kmen1=menu('MAKE OUTPUT FILES (choose type)',...
                '.RWL',...
                '.RW, .EWW or .LWW',...
                'QUIT');
            if kmen1==1; % .rwl files
                kmen = menu('Choose type of .rwl file',...
                    'TOTAL RING WIDTH                ',...
                    'EARLYWOOD WIDTH          ',...
                    'LATEWOOD WIDTH  ',...
                    'OTHER                       ');
                if kmen==1;
                    ftype='rw';
                    rwc2rwl(XT,XEL,pathout,ftype,rwlset);
                elseif kmen==2; % .rwl files of early widtH
                    ftype='ew';
                    rwc2rwl(XT,XEL,pathout,ftype,rwlset);
                    %uiwait(msgbox('NOT IMPLEMENTED YET','Message','modal'));
                elseif kmen==3; % .rwl files of late widt
                    ftype='lw';
                    rwc2rwl(XT,XEL,pathout,ftype,rwlset);
                    %uiwait(msgbox('NOT IMPLEMENTED YET','Message','modal'));
                elseif kmen==4 % .rwl file of total width - measured only
                    kmen2=menu('.RWL FILE OF TOTAL WIDTH',...
                        'FROM TOTAL-RING MEASUREMENTS ONLY)',...
                        'FROM PARTIAL-RING MEASUREMENTS ONLY');
                    if kmen2==1;% from measure only;
                        ftype='rwm';
                        rwc2rwl(XT,XEL,pathout,ftype,rwlset); 
                    elseif kmen2==2; % from computed only
                        ftype='rwc';
                        rwc2rwl(XT,XEL,pathout,ftype,rwlset);
                        %uiwait(msgbox('NOT IMPLEMENTED YET','Message','modal'));
                        
                    end;
                    clear kmen2;
                end; % if kmen==1;
                clear kmen;
            elseif kmen1==2; % .rw, .eww or .lww files
                kmen = menu('CHOOSE OUTPUT',...
                    '.RW FILES OF TOTAL RING WIDTH          ',...
                    '.EWW & .LWW FILES OF EARLYWOOD AND LATEWOOD WIDTH          ',...
                    'OTHER                       ');
                if kmen==1; % .rw files
                    ftype='rw'; % total width from XEL and XT if available
                    rwc2rw(XT,XEL,pathout,ftype);
                elseif kmen==2; % .eww & .lww files
                    ftype='ew/lw';
                    rwc2rw(XT,XEL,pathout,ftype);
                elseif kmen==3; % other
                    kmen2=menu('.RW FILE OF TOTAL WIDTH FROM:',...
                        'FROM TOTAL-RING MEASUREMENTS ONLY)',...
                        'FROM PARTIAL-RING MEASUREMENTS ONLY');
                end;
            elseif kmen1==3;
                kwh1=0;
            end; % if kmen1
        end; % while kwh1
    end; % if kedit
    switch kedit1;
    case{1,2,3,4}; % updatind needed after editing
         
        % Update span
        
        if kedit1 ~=4; % if not just modifying an rwlset
            eval([vbl ' = subfcn07(' vbl ');']);
            
            % Sort the revised data so that series alphabetically ordered
            eval([vbl ' = subfcn02(' vbl ');']);
            
            % Update X.summary (string)
            if strcmp(editmode,'Total');
                eval([vbl ' = subfcn03(' vbl ',1);']); % the "1" means total ring width
            elseif strcmp(editmode,'Partial');
                eval([vbl ' = subfcn03(' vbl ',2);']);
            else; 
                error ('Invalid editmode');
            end;
            
            % Update string info -- sinfo
            sinfo=subfcn05(XT,XEL);
        end;
        
        % Save results -- call subfcn04
        subfcn04(filemode,pf1,XT,XEL,vlist,rwlset,sinfo);
        
    otherwise; % wrote output files, no updating needed
    end;
        
case 'Measure';  % MEASUREMENT MODE
    
    kmen = menu('Choose','Start new site',...
        'Operate on existing site');
    switch kmen;
    case 1; 
        filemode='New';
        [XT,XEL,sinfo,rwlset,vlist,pf1]=subfcn08(units);
    case 2;
        filemode='Existing';
        %-- Prompt for  name of existing rw storage file
        [file1,path1]=uigetfile('*.mat','Input file storing ring-width data');
        pf1=[path1,file1];
        if ~exist(pf1,'file'); % check that file exists
            error([pf1 ' does not exist']);
        else; % File exists; check that it has the required data
            eval(['s=load(''' pf1 ''');']);
            if ~all([isfield(s,'vlist')  isfield(s,'XEL')  isfield(s,'XT') isfield(s,'sinfo')  ]);
                error([pf1 ' does not contain XEL, XT vlist and sinfo']);
            end
            clear s;
            % Load the ring-width data
            eval(['load ' pf1 '  XT XEL sinfo vlist rwlset;']); 
        end; % if ~exist(pf1,'file');
        
    otherwise;
    end;
    
    
    %--- PARTIAL OR TOTAL RING WIDTH 
    
    kmen = menu('Choose Ring-Width Type',...
        'Total Ring Width',...
        'Earlywood/Latewood Width');
    switch kmen;
    case 1; 
        widthmode='Total';
        vbl='XT';
        jcol = [2]; % column of XT.data holding time series
    case 2;
        widthmode='EW/LW';
        vbl='XEL';
        jcol = [2 3 4]; %% columns of XEL.data holding time series
    otherwise;
    end;
    
    
    %--- WHICH SERIES? Get the series id and its sequence number
    
    if strcmp(filemode,'New');
        kseries='New'; % new series, first of a new file
        jser=1; % first series
        nameser = subfcn01; % Call subfcn01.m to  name series
    else; % filemode not 'New'
        kmen = menu('Choose',...
            'Start measuring a new series',...
            'Continue measuring an existing series');
        if kmen==1; % start measuring a new series in an existing storage file
            kseries = 'New';
            eval(['idany = ' vbl '.id;']);
            
            if isempty(idany{1});
                jser=1;
            else;
                jhi=length(idany);
                jser = jhi+1;
            end;
            yrsagain=[];
            
            kwh2=1; % while loop for naming series
            while kwh2==1;
                nameser = subfcn01; % Call subfcn01.m to  name series
                % Check that this "new" series does not have same name as an existing series
                if any(strcmp(idany,nameser));
                    uiwait(msgbox([nameser ' is already a series in the measurement file'],'Invalid id for new series','modal'));
                else;
                    kwh2=0;
                end;
            end; % while kwh2==1;
            clear kwh2;
            
            
        else;  % resume measuring an existing series -- or remeasure parts
            kseries='Existing';
            eval(['xid = ' vbl '.id;']);
            eval(['yrsa = ' vbl '.remeasure;']);
            
            kmen = menu('Choose series to measure',xid);
            if ~isempty(yrsa{kmen});
                strwarn =num2str(yrsa{kmen});
                beep; beep; beep;
                uiwait(msgbox(strwarn,'Remeasure these years first','modal'));
            end;
                
            jser=kmen;
            nameser=xid{jser};
            clear kmen xid;
        end
    end; % if strcmp(filemode,'New');
    
    
    
    %--- SET SOME YEAR INFORMATION FOR THE SERIES
    
    % SPECIFY yrorigin, yrstart, YRSTART, AND INITIALIZE xin and yrin.
    % yrorigin is the first year of the STORED width series.  yrstart is the year you want to begin 
    % measuring at.  xin is a time series vector (if vbl is XT) or matrix (if vbl is XEL) of current
    % measurements to be passed to the measuring function. xin might be filled with NaN if no measurements yet 
    % exist. xin will be NaN filled to a length maxlen years before passing it to the measuring function.
    
    % Set yrorigin and yrstart
    
    if strcmp(filemode,'New') | strcmp(kseries,'New');
        % Prompt for first year of measurement
        prompt={'Enter year of first value to be measured:'};
        def={'1'};
        dlgTitle='Input Year';
        lineNo=1;
        answer=inputdlg(prompt,dlgTitle,lineNo,def);
        yrstart = str2num(answer{1});
        yrorigin=yrstart-nyrlead;  
        yrsagain=[];
        
    else; % the series is pre-existing
        eval(['xdata = ' vbl '.data;']);
        eval(['remeas = ' vbl '.remeasure;']);
        xdata=xdata{jser};
        yrsagain=remeas{jser};
       
        yrfirst= min(xdata(:,1)); % first year of existing measurements tsm
        yrorigin=yrfirst-nyrlead;
        yrlast=  max(xdata(:,1));  % last ...
        
        if isempty(yrsagain);
            yrstart = yrlast+1;
        else;
            yrstart=yrsagain(1);
        end;
        
        strmeas1=['Series ' nameser ', ' num2str(yrfirst) '-' num2str(yrlast)];
        strmeas2=['Start measuring at year '];
        yrdef=yrstart;
        prompt={strmeas2};
        def={num2str(yrdef)};
        dlgTitle=strmeas1;
        lineNo=1;
        answer=inputdlg(prompt,dlgTitle,lineNo,def);
        yrstart=str2num(answer{1})
                
    end;
    clear prompt def dlgTitle lineNo answer; 
    
    
    
    %---- INITIALIZE THE DATA ; Initialize xin and yrin
    
    if strcmp(widthmode,'Total');
        xin = repmat(NaN,maxlen,1);
    else;
        xin = repmat(NaN,maxlen,3);
    end;
    
    yrsp = yrorigin+maxlen-1;
    yrin=(yrorigin:yrsp)';
    
    % Postion the existing measurements in xin
    if strcmp(filemode,'New') | strcmp(kseries,'New');
        % no action needed, xin is all NaN
    else;
        igo = yrfirst-yrorigin+1;
        isp = yrlast-yrorigin+1;
        xdata(:,1)=[];
        xin(igo:isp,:)=xdata;
        if size(xdata,1)>maxlen-2*nyrlead;
            strmess1={['Possible problem with too long a series'],...
                    ['Series length is ' int2str(size(xdata,1))],...
                    ['while hard coded settings allow for ' int2str(maxlen) 'yr,'],...
                    ['which includes ' int2str(nyrlead) ' leading years.'],...
                    'You might run out of room in adding additional measurements.',...
                    'Solution: look into changing hard coded defaults for ',...
                    'maxlen and nyrlead at start of this function'};
           uiwait(msgbox(strmess1,'Warning','modal'));
       end;
    end;
    
    
    
    % -- MEASURE
    
    if strcmp(widthmode,'Total');
        [xout,yrout,yrsagain]=meastot(xin,yrin,yrstart,yrsagain); 
    else;
        [xout,yrout,yrsagain]=measpart(xin,yrin,yrstart,yrsagain);
    end;
    
    
    ddate=date; %  compute current date
    
    % Store measurements, span, etc
    
    eval([vbl '.id{' int2str(jser) '} =nameser;']);
    eval([vbl '.data{' int2str(jser) '} = [yrout xout];']);
    eval([vbl '.who{' int2str(jser) '} = person;']);
    eval([vbl '.when{' int2str(jser) '} = ddate;']);
    eval([vbl '.remeasure{' int2str(jser) '} = yrsagain;']);
    if strcmp(widthmode,'EW/LW');
        % Compute spans for Ew and Lw; store spans
        iE = find(~isnan(xout(:,1)));
        iL = find(~isnan(xout(:,2)));
        iC = find(~isnan(xout(:,3)));
        eval([vbl '.span{' int2str(jser) '} = [yrout(min(iE))  yrout(max(iE)); yrout(min(iL)) yrout(max(iL)) ; yrout(min(iC))  yrout(max(iC)) ];']);
    elseif strcmp(widthmode,'Total');
        eval([vbl '.span{' int2str(jser) '} = [min(yrout) max(yrout)];']);
    else;
        error('Invalid widthmode');
    end;
    
    
    % Sort the revised data so that series alphabetically ordered
    eval([vbl ' = subfcn02(' vbl ');']);
    
    % Update X.summary (string)
    if strcmp(widthmode,'Total');
        eval([vbl ' = subfcn03(' vbl ',1);']); % the "1" means total ring width
    else;
        eval([vbl ' = subfcn03(' vbl ',2);']);
    end;
    
    % Update string info -- sinfo
    sinfo=subfcn05(XT,XEL);
    
    % Save results -- call subfcn04
    subfcn04(filemode,pf1,XT,XEL,vlist,rwlset,sinfo);
    
   
otherwise;
end; %  switch runmode;





% ******  SUBFUCTION 1  -- getting and checking the series name

function nameser = subfcn01;

prompt={'Enter the series ID:'};
def={'XXX01A'};
dlgTitle='Series ID (site-code + tree# + core=letter (e.g., PDF01A)';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
d1=answer{1};
d2=deblank(d1);
d3=fliplr(deblank(fliplr(d2)));
nameser=upper(d3);
clear d1 d2 d3 def dlgTitle lineNo answer;

% Length of name
nlen = length(nameser);
if nlen>8;
    error(['Series ID ' nameser ' has more than 8 chars']);
elseif nlen<5;;
    error(['Series ID ' nameser ' needs at least 5 chars for a site code, tree no and core letter']);
    
end;


% First 3 chars must be letters
d1 = nameser(1:3);
if ~all(isletter(d1));
    error([nameser ' first 3 chars must be letters']);
end;
sitecode=d1;
clear d1;

% Check that have a core letter, following a tree number
d1=nameser;
d1(1:3) = []; % strip off site code
i1 = find(isletter(d1)); % find remaining letters
if isempty(i1);
    error([nameser ' needs a tree number following site code']);
else;
    treeno=d1(1:i1(1)-1);
    if length(treeno)==1;
        treeno=['0' treeno];
    elseif length(treeno)>3;
        error([nameser ' tree number should be max of 999 and max of 3 digits']);
    end;
end;

% Check core letter and for trailing core-segment number
d1 = d1(i1(1):length(d1));
i1=find(isletter(d1));
if length(d1)==1;
    coreletter=d1;
    nameser=[sitecode treeno coreletter];
else;
    % Strip the letter off, and should have only a 1-digit number
    coreletter=d1(i1(1));
    d1(1)=[];
    if length(d1)>1;
        error([nameser ' has a core-segment number after the core letter, but this should be length 1']);
    else;
        if isletter(d1); 
            error(['Last char of ' nameser ' should be a number for this case']);
        end;
        segno = d1;
        nameser=[sitecode treeno coreletter segno];
    end;
end;


% SUBFUNCTION 2 -- sort the cell data

function [X]=subfcn02(X)
% Sort cell storage alphabetically

% id, data , when, who, span
if length(X.id)>1;
    [X.id,j]=sort(X.id);
    X.data=X.data(j);
    X.when=X.when(j);
    X.who=X.who(j);
    X.span=X.span(j);
    X.remeasure=X.remeasure(j);
else;
    % No action needed to sort one
end;


%--- SUBFUNCTION 3 --  update the summary
function X=subfcn03(X,k)
% update the summary

n= length(X.id); % number of series measured

if k==1; % total rw
    str1='Total Ring Width';
    str2=' ID         Total RW';
elseif k==2; % EWW/LWW
    str1='EWW/LWW';
    str2=['ID            EWW           LWW           EWW+LWW'];
end;

s1=[int2str(n) ' series  have been measured for ' str1];
s1=char(s1,'Spans of measurements are as follows');
s1=char(s1,str2);


for j = 1:n;
    spann=X.span{j};
    if k==1;;
        s2=[X.id{j} '      '  int2str(spann(1,:))];
    elseif k==2;
        s2 = [X.id{j} '     '   int2str(spann(1,:)) '     ' int2str(spann(2,:)) '     '  int2str(spann(3,:))];
    end;
    s1=char(s1,s2);
end;

X.summary=s1;


%-------  SUBFUNCTION 4  UPDATING FILES

function subfcn04(filemode,pf1,XT,XEL,vlist,rwlset,sinfo);
set1=' XT XEL sinfo rwlset vlist';
switch filemode;
case 'Existing';
    kmen=menu('Choose',...
        ['Save revised data in ' pf1],...
        'Save revised data in a new file',...
        ['Abort -- original ' pf1 ' will be unchanged']);
    if kmen==1;
        eval(['save ' pf1  set1 ' -append;']);
    elseif kmen==2;
        [file2,path2]=uiputfile('*.mat','Output storage file');
        pf2=[path2 file2];
        eval(['save ' pf2  set1 ';']);
    elseif kmen==3;
        disp('ABORTED ');
        return;
    end;
case 'New';
    kmen=menu('Choose',...
        ['Save the data in new storage file ' pf1],...
        ['Abort ']);
    if kmen==1;
        eval(['save ' pf1  set1 ';']);
    elseif kmen==2;
        disp('ABORTED ');
        return;
    end;
   
otherwise; % switch filemode -- saving file
end;% switch filemode -- saving file

%---- SUBFUNCTION 5  UPDATING THE STRING SUMMARY

function sinfo=subfcn05(XT,XEL)
% summary:  summary information on total and partial ring-width series
% sinfo=summary(sinfo,X,XE);
% Last revised 10-03-01
%
% Summary text info on measurements stored
%
%*** INPUT
%
% XT, XEL  structures of total and partial ring width (see rwmeas.m)
%
%
%*** OUTPUT
% 
% sinfo -- character array with following cols
%   id -- with * if remeasure flagged for any variable
%   total width
%       start and end year for composite
%       start and end year for measured only
%       start and end year for computed only
%   EWW
%       start and end year
%   LWW 
%       start and end year
%
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED
% uniqcell
%
%*** TOOLBOXES NEEDED-- NONE
%
%*** NOTES
%


% Store the structure data
id1=XT.id; % series ids
id2=XEL.id;
d1=XT.data; % series time series, vector or matrix
d2=XEL.data;
r1=XT.remeasure; % remeasure years
r2=XEL.remeasure;

idset = uniqcell(id1,id2); % get cell array of series ids, no duplicates, from XT and XEL
nset = size(idset,2); % number of ids

% Initialize

D=repmat(NaN,nset,10);

for n =1:nset;
    idthis=idset(n);
    
    % Any rw series?
    j1 = strcmpi(idthis,id1);
    if ~all(j1==0); % measured total width exists for this series
        X1=d1{j1}; % total width data
        yrx1 = X1(:,1);
        x1=round(100*X1(:,2));
        % Strip trailing and leading NaN
        [x1,yrx1]=subfcn06(x1,yrx1);
        D(n,1)=yrx1(1);
        D(n,2)=yrx1(end);
        D(n,3)=yrx1(1);
        D(n,4)=yrx1(end);
    else;
        % No action; D entries remain NaN
    end;
    
    % Any eww/lwwseries?
    j2 = strcmpi(idthis,id2);
    if ~all(j2==0); % measured total width exists for this series
        X2=d2{j2}; % total width data
        yrx2 = X2(:,1);
        xearly=round(100*X2(:,2));
        xlate=round(100*X2(:,3));
        xcomp= round(100*X2(:,4));
        % Strip trailing and leading NaN
        [xearly,yrxearly]=subfcn06(xearly,yrx2);
        [xlate,yrxlate]=subfcn06(xlate,yrx2);
        [xcomp,yrxcomp]=subfcn06(xcomp,yrx2);
        
        % merged total
        D(n,1)=min([yrxcomp(1) D(n,1)]);;
        D(n,2)=max([yrxcomp(end) D(n,2)]);
        
        % Early
        D(n,7)=yrxearly(1);
        D(n,8)=yrxearly(end);
        
        % Late
        D(n,9)=yrxlate(1);
        D(n,10)=yrxlate(end);
        
        % Computed total
        D(n,5)=yrxcomp(1);
        D(n,6)=yrxcomp(end);
        
    else;
        % No action; D entries remain NaN
    end;
end;

D=num2str(D,'%8.0f');
D=[char(idset) repmat(blanks(2),nset,1) D];

hd1='                      TOTAL WIDTH                             PARTIAL WIDTH         ';
hd2='        --------------------------------------------   -----------------------------';
hd3='           MERGED        MEASURED      COMPUTED          EARLYWOOD        LATEWOOD  ';
hd4='ID      --------------------------------------------    ------------    ------------';

sinfo = char(hd1,hd2,hd3,hd4,D);
        
%---- SUBFUNCTION 06    
function [x,yrx]=subfcn06(x,yrx);
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


%--- SUBFUNCTION 07  UPDATE .SPAN
function X=subfcn07(X);
nser = length(X.id);
for n =1:nser;
    x=X.data{n};
    yrx=x(:,1);
    if size(x,2)==4; % partial
        iE = find(~isnan(x(:,2)));
        iL = find(~isnan(x(:,3)));
        iC = find(~isnan(x(:,4)));
        X.span{n}=[yrx(min(iE))  yrx(max(iE)); yrx(min(iL)) yrx(max(iL)) ; yrx(min(iC))  yrx(max(iC)) ];
    elseif size(x,2)==2;  % total 
        X.span{n}= [min(yrx) max(yrx)];
    else;
        error ('Invalid col size for x');
    end;
end; 



%---- SUBFUNCTION 08   INITIALIZE STORAGE FOR A NEW ,MAT FILE

function [XT,XEL,sinfo,rwlset,vlist,pf1]=subfcn08(units);
% Initialize variables for a new site
%-- Prompt for  name of new storage file
[file1,path1]=uiputfile('*.mat','New input file to store ring widths in ');
pf1=[path1,file1];
if exist(pf1,'file'); % check that file exists
    error([pf1 ' already exists']);
else; % File exists; check that it has the required data
end; % if ~exist(pf1,'file');

% Initialize variables
XT.id{1}=[]; % sample id -- same for total, EWW or LWW, but empty if does not exist 
XEL.id{1}=[];
XT.data{1}=[];
XEL.data{1}=[];
XT.span{1}=[];
XEL.span{1}=[];
XT.who{1}=[];
XEL.who{1}=[];
XT.when{1}=[];
XEL.when{1}=[];
XT.remeasure{1}=[];
XEL.remeasure{1}=[];
XT.units=units;
XEL.units=units;
XT.summary=[];
XEL.summary=[];
rwlset.name{1}=[];
rwlset.describe{1}=[];
rwlset.when{1}=[];
rwlset.trimall{1}=[];
rwlset.idnames{1}=[];
rwset.trimeach{1}=[];
sinfo=[]; % string information on individual series


vlist = 'Contents of .mat file';
vlist =char(vlist,'');
vlist=char(vlist,'XT,XEL - structure variables of total ring width and earlywood/latewood width:');
vlist =char(vlist,'   XL.units, XEL.units -- units of measurement');
vlist =char(vlist,'   XL.summary, XEL.summary',...
    '   -- series id,  first and last year of data (string matrix)');
vlist =char(vlist,'   NEXT VARIABLES ARE CELLS, ONE ELEMENT PER SERIES, for XT, XEL');
vlist =char(vlist,'     {id -- name of the series (e.g., PDF02A)}');
vlist =char(vlist,'     {data -- the measurements}');
vlist =char(vlist,'     {span -- first and last year of measurements}');
vlist =char(vlist,'     {who -- the initials of measurer}');
vlist =char(vlist,'     {when -- date measurements last modified}');
vlist =char(vlist,'     {remeasure -- years flagged to remeasure}');
vlist =char(vlist,'sinfo - string summary of span of all types of data for unique series names');
vlist =char(vlist,'rwlset -- cell arrays of ids defining sets of series to be used in making rwl files');
vlist =char(vlist,'     {name -- string name (max of 20 chars of rwl set }');
vlist =char(vlist,'     {when -- date the rwlset created}');
vlist =char(vlist,'     {describe -- description of the rwl set (the desciption is a cell array of strings}');
vlist =char(vlist,'     {trimall -- restricted time coverage for the rwlset (first and last year); [] if no restriction}');
vlist =char(vlist,'     {idnames -- the id names of series to be included in the rwl set}');
vlist =char(vlist,'     {trimeach -- restricted time coverage for individual series in the rwlset; [] if no restriction}');

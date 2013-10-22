function rwlset=rwlspecs(rwlset,XT,XEL);
% rwlspecs:  make or modify list of ids for use in a .rwl file
% rwlset=rwlspecs(rwlset);
% Last modified 2006-10-11
%
% Make or modify list of ids for use in a .rwl file.  You want to automate the process of
% selecting series for an .rwl file so that you can regenerate the file at any time. This
% function makes a list you can later recall in writing .rwl files with rwmeas.m
%
%*** INPUT 
%
% rwlset.  structure of .rwl file sets, with fields:
%   .name {j}s    short name or code of rwlset (e.g., padwt1}
%   .describe{j}   cell array of text description of the jth .rwl set
%   .trimall{j}    (1 x 2)i specified time coverage of data in the rwlset (all series truncated to this)
%   .idnames{j}  cell array of id names in jth set
%   .trimeach{j}  cell array of specified start and end years of spans for indiv series in the rwlset
% XT structure variable of total width
% XEL structure variable of partial widths
%
%*** OUTPUT
%
% rwlset.  revised structure
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% .describe{}  rwlset.describe{1} holds the description for the first id set, or is [] if no
%   sets defined.  rwlset.describe{1} is itself a cell array of strings
% .idnames{}  rwlset.idnames{1} holds a cell array of series ids for series to be included in the
%   rwl set .  For example rwlset.idnames{1} ={'pdf04a','pdf05a'}
% Origin.  rwlspecs.m first written as function to be called by rwmeas.m in measurement program
%
% Options
% Possible options are
%   1) create a new rwlset, regular
%   2) remove
%   3) modify
%   4) create a latewood rwlset interactively with viewing of plots
%


% Except for latewood specials, will only need id1 and id2 
id1=XT.id; %cell array of total width (XT) names currently in the storage file
id2=XEL.id; %cell array of EWW (also LWW) names currently in the storage file
clear XT;


% Numbers of series currently stored
num1=length(id1); % number of XT ids
if num1==1 & isempty(id1{1});
    num1=0;
end;
num2=length(id2); % number of XEL ids
if num2==1 & isempty(id2{1});
    num2=0;
end;

% -- HOW MANY RWLSETS NOW EXIST

if isempty(rwlset.name{1});
    nsets=0;
else;
    nsets=length(rwlset.name);
end;
nextset=nsets+1; % next set, if creating new set



% -- OPTION FOR CREATING, DELETING, MODIFYING RWL SET

kmen=menu('Choose Option',...
    'Create a new rwlset, regular',...
    'Remove an rwlset',...
    'Modify an existing rwlset',...
    'Create a latewood rwlset interactively');
if kmen==1; % 'Create a new rwlset'
    kmode='New';
elseif kmen==2; % remove an rwlset
    kmode='Remove';
elseif kmen==3;
    kmode='Modify'; % modify an existing rwl set by marking series to be added or deleted
elseif kmen==4;
    kmode='LW'; %  'Create a latewood rwlset interactively
end;
clear kmenl; 


%---- IF LATEWOOD SPECIAL, HANDLE THAT NOW AND THEN RETURN TO CALLING FUNCTION

if strcmp(kmode,'LW');
     rwlset=rwlspecsLW(rwlset,XEL);
     return;
end



% -- WHAT TYPE OF SERIES WILL THE .RWL BE USED FOR -- THIS DETERMINES THE SET OF CANDIDATE SERIES
%  TO ICNLUDE

kmen = menu('Choose type of target .rwl series',...
    'TOTAL RING WIDTH                ',...
    'EARLYWOOD WIDTH OR LATEWOOD WIDTH ',...
    'OTHER                       ');
if kmen==1;
    ftype='rw';
    % Build default menus
    idset=uniqcell(id1,id2);
    numset=length(idset);
    yesx = cellstr(repmat('-Y',numset,1));
    nox  = cellstr(repmat('-N',numset,1));
    grpyx = cellstr([char(idset) char(yesx)]);
    grpnx = cellstr([char(idset) char(nox)]);
        
elseif kmen==2; % .rwl files of earlywood and latewood width
    ftype='ew';
    
    % Build default menus
    yesx = cellstr(repmat('-Y',num2,1));
    nox = cellstr(repmat('-N',num2,1));
    idset=id2;
    numset=length(idset);
    grpyx = cellstr([char(idset) char(yesx)]);
    grpnx = cellstr([char(idset) char(nox)]);
    
elseif kmen==3 % special case series of total width
    kmen2=menu('SPECIAL .RWL FILE OF TOTAL WIDTH',...
        'FROM TOTAL-RING MEASUREMENTS ONLY)',...
        'FROM PARTIAL-RING MEASUREMENTS ONLY');
    if kmen2==1;% from measure only;
        ftype='rwm';
        % Build default menus
        idset=id1;
        numset=length(idset);
        yesx = cellstr(repmat('-Y',numset,1));
        nox  = cellstr(repmat('-N',numset,1));
        grpyx = cellstr([char(idset) char(yesx)]);
        grpnx= cellstr([char(idset) char(nox)]);
        
    elseif kmen2==2; % from computed only
        ftype='rwc';
        % Build default menus
        idset=id2;
        numset=length(idset);
        yesx = cellstr(repmat('-Y',numset,1));
        nox = cellstr(repmat('-N',numset,1));
        grpyx = cellstr([char(idset) char(yesx)]);
        grpnx = cellstr([char(idset) char(nox)]);
        
    end;
    clear kmen2;
end; % if kmen==1;


% -- GO THRU MODES

nallow=[1 numset]; % must select at least 1 series, no more than numset 

if strcmp(kmode,'New'); % -- MAKE NEW RWLSET
    
    kmen=menu('Choose one',...
        'Start with full plate',...
        'Start with empty plate');
    if kmen==1;
        kplate='Full';
        strmenu='Toggle Y to N to remove a series';
        C1=idset; % cell array of names
        C2=yesx;
        Lin = logical(ones(numset,1));
        [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu); % C3 is cell array of selected series
    else; 
        kplate='Empty';
        strmenu='Toggle N to Y to add a series';       
        C1=idset; % cell array of names
        C2=nox;
        Lin = logical(zeros(numset,1));
        [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu); % C3 is cell array of selected series
    end;
    rwlset.idnames{nextset}=C3; % store the ids for this rwlset
    
    % Initially, specify no restriction on the time coverage of the rwlset
    rwlset.trimall{nextset}=[];
    
    % Initially set no restrictions on the individual spans to be included in the rwlset
    rwlset.trimeach{nextset}=cell(1,size(C3,2));

    %-- Name the new rwlset
    prompt={'Enter short name for rwl set','Enter description:'};
    def={'DefaultSet','Default set of series to send in .rwl file to ITRDB'};
    dlgTitle='Name and description for this rwlset';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    rwlset.name{nextset}=answer{1};
    rwlset.describe{nextset}=answer{2};
    
    % Store the creation date
    rwlset.when{nextset}=date;
    
      
elseif strcmp(kmode,'Remove'); % remove an rwlset
    
    if nsets==0;
         uiwait(msgbox('No rwl sets to remove','Message','modal'));
         return
    else;
        C1=(rwlset.name)';
        C2=cellstr(repmat('-N',nsets,1));
        if isempty(C1{end});
            C1(end)=[];
            C2(end)=[];
            nsets = length(C1);
            if isempty(C1);
                error('No rwl sets');
            end
        end
            
        Lin = logical(zeros(nsets,1));
        nallow=[1 nsets];
        mess2= ' Toggle Y to N in following menu to delete any rwlset(s)';
        uiwait(msgbox(mess2,'Message','modal'));
        strmenu='Toggle Y to N to remove any rwlset(s)'
        [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu);
        rwlset.name(find(Lout))=[];
        rwlset.describe(find(Lout))=[];
        rwlset.idnames(find(Lout))=[];
        rwlset.when(find(Lout))=[];
        rwlset.trimall(find(Lout))=[];
        rwlset.trimeach(find(Lout))=[];
        if all(Lout);
            rwlset.name{1}=[];
            rwlset.describe{1}=[];
            rwlset.idnames{1}=[];
            rwlset.when{1}=[];
            rwlset.trimall{1}=[];
            rwlset.trimeach{1}=[];
        else;
            
        end;
    end;
      
    
elseif strcmp(kmode,'Modify');
    
    if nsets==0; 
        error('Cannot modify rwlset  because no rwlsets yet created');
    end;
        
    %--CHOOSE THE RWLSET TO MODIFY
    
    nameset=rwlset.name;
    noxx = cellstr(repmat('-N',nsets,1));
    grpnxx = cellstr([char(nameset) char(noxx)]);
    strmenu='Toggle N to Y to select rwlset to modify';
    C1=nameset; % cell array of names of rwlsets
    C2=noxx;
    Lin = logical(zeros(nsets,1));
    nallow=[1 1];
    [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu); % C3 is cell array of selected series
    
    %--- PULL THE NAMES FOR THE SET TO BE MODIFIED
    iset =find(Lout); % the selected rwlset
    namethis = nameset{iset}; %  name of the rwlset to be modified
    idsthis = rwlset.idnames{iset}; %  cell of the names of the series currently in the rwlset
    trimthis=rwlset.trimeach{iset};
    
    %-- BUILD MENU FOR MODIFYING SERIES
    
    L=cellmch1(idsthis,idset);
    islot = cellmch2(idsthis,idset);
    if all(isnan(islot));
        error('No matches of rwlset ids with available ids');
    else;
        nfound=sum(~isnan(islot));
    end;
    if nfound ~= sum(L);
        error('Mismatch of number of non-NaN elements of islot with sum of 1s in L');
    end;
    trimworld = cell(1,length(idset));
    if ~isempty(trimthis); % if any of the rwlset series had a trimmed period
        trimworld(L)=trimthis(islot(~isnan(islot))); % cell of trim periods (if any) associated with each id in idset
    else; % Associate [] as starting trimmed period of each series in rwlset, since none had a specified trim period
        % no action needed.  trimworld(L) already []
    end;
       
    ithis =find(L);
    numset=length(idset);
    nox  = cellstr(repmat('-N',numset,1));
    nox(ithis)={'-Y'};
    grpnx = cellstr([char(idset) char(nox)]);
    Lin=L;
    
    %-- CHOOSE HOW TO MODIFY
    
    kmen = menu(['Choose how to modify rwlset ' namethis ],...
        'Add or remove a series',...
        'Restrict the period for some individual series in the rwl files',...
        'Restrict the period for the entire rwl file');
    if kmen==1; %  Add or remove a series
        strmenu='Toggle Y to N for series to be removed from rwlset';
        nallow=[1 numset];
        C1=idset;
        C2=nox;
        [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu); % C3 is cell array of selected series
        ipicked=find(Lout);
        % Update the cell of series in the rwlset
        rwlset.idnames{iset}=C3;
        rwlset.trimeach{iset}=trimworld(ipicked);
    elseif kmen==2; % Restrict the period for some individual series in the rwl files
        mess1 ={'You will enter a menu to select series for which you want to truncate, or restrict',...
                'the period to be written to the .rwl file.  In the series-selection menu, series',...
                'whos period is restricted are followed by the start and end year of the',...
                'restricted period (e.g., [1882 1989]).  If there is no restriction, no years follow',...
                'the series id in the menu'};
        uiwait(msgbox(mess1,'Message','modal'));
        kwh1=1;
        
        
        % Build a char matrix of the restricted period (or [-inf inf]) for series
        navail=size(idsthis,2); % number of series in this rwset
        % Make a matrix of the restricted periods
        P = repmat([-inf inf],navail,1);
        for n=1:navail;
            p = rwlset.trimeach{iset}{n};
            if ~isempty(p);
                P(n,:)=p;
            else;
            end;
        end;
        P0=P;
        Pstr=num2str(P0);
        C1=cellstr([char(idsthis) repmat(blanks(1),navail,1)   Pstr  repmat(blanks(1),navail,1)]);
        
        while kwh1;
            
            
            kmen1=menu('Choose one','Change restriction on output for another series',...
                'Finished changing restrictions');
            if kmen1==1; % Change restriction on output for another series
                strmenu='Toggle N to Y for the series you want to reset the output period for';
                nallow=[1 1]; % operate on series at a time
                noox  = cellstr(repmat('-N',navail,1));
                C2=noox;
                Lin=logical(zeros(navail,1));
                [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu); % C3 is cell array of selected series
                ipicked=find(Lout);
                
                % Set the restriction
                prompt={['Enter the period (' C3{1} '):']};
                
                kwh2=1;
                while kwh2;
                   
                    def={num2str(P0(ipicked,:))};
                    dlgTitle=['Restrict rwl dat for  ' rwlset.idnames{iset}{ipicked} ' to'];
                    lineNo=1;
                    answer=inputdlg(prompt,dlgTitle,lineNo,def);
                    p0=str2num(answer{1});
                    xdiff=diff(p0);
                    if ~isnumeric(xdiff) |  length(xdiff)~=1 | isnan(xdiff) | xdiff==-inf | xdiff<=0;
                        mess2=[num2str(p0) ' invalid as start and end of a period. Do it again'];
                        uiwait(msgbox(mess2,'Message','modal'));
                    else;
                        kwh2=0;
                    end;
                end;
                P0 (ipicked,:)= p0;
                Pstr=num2str(P0);
                C1=cellstr([char(idsthis) repmat(blanks(1),navail,1)   Pstr  repmat(blanks(1),navail,1)]);
            elseif kmen1==2; % Finished changing restrictions
                kwh1=0;
            end; % menu('Choose one','Change restriction...
        end;  % while kwh1;
        
        % Store the revised restricted periods
        for n=1:navail;
            p=P0(n,:);
            if p(1)==-inf & p(2)==inf;
                p=[];
            else;
            end;
            rwlset.trimeach{iset}{n}=p;
        end;
    elseif kmen==3; % Restrict the period for the entire rwl file
        % Recall the C3 is the name of the selected rwlset, and that Lout is logical pointer to it
        
        kwh3=1;
        g =     rwlset.trimall{iset};
        if isempty(g);
            g0=[-inf inf];
        else;
            g0=g;
        end;
        while kwh3;
            ithis =find(Lout);
            
            prompt={'Enter first and last year (-inf and inf also acceptable):'};
            def={num2str(g0(1,:))};
            dlgTitle=['Restrict the output rwlset ' rwlset.name{iset} ' to period:'];
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            g0=str2num(answer{1});
            xdiff=diff(g0);
            if ~isnumeric(xdiff) |  length(xdiff)~=1 | isnan(xdiff) | xdiff==-inf | xdiff<=0;
                mess2=[num2str(g0) ' invalid as start and end of a period. Do it again'];
                uiwait(msgbox(mess2,'Message','modal'));
            else;
                kwh3=0;
            end;
        end;
        % Store revised period
        g=g0;
        if g(1)==-inf & g(2)==inf;
            g=[];
        else;
        end;
        rwlset.trimall{iset}=g;
    end; % menu(['Choose how to modify rwlset '
            
        
    %-- CODE TO CORROBORATE CHANGE
    
    %-- CODE TO UPDATE .IDNAMES ,.WHEN, .COVER, .SPAN AS APPLICABLE
    
    % Update the name and description of the rwlset
    prompt={'Enter short name for rwl set','Enter description:'};
    def={rwlset.name{iset},rwlset.describe{iset}};
    dlgTitle='Edit the name and description for this rwlset if you want';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    rwlset.name{iset}=answer{1};
    rwlset.describe{iset}=answer{2};
    
    % Store the revision date
    rwlset.when{iset}=date;
else;
    error('Invalid kmode');
end;


%--- CLEAN UP POSSIBLE PROBLEM WHEN HAVE AN EMPTY RWLSET AFTER A REAL ONE
if (length(rwlset.name)>1 & isempty(rwlset.name{end}));
    rwlset.name(end)=[];
    rwlset.describe(end)=[];
    rwlset.when(end)=[];
    rwlset.idnames(end)=[];
    rwlset.trimall(end)=[];
    rwlset.trimeach(end)=[];
end
       





    

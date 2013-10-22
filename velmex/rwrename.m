function [XT,XEL,rwlset,kboth]=rwrename(XT,XEL,rwlset,editmode)
% rwrename:  rename, or assign new ids to, one or more ring-width series
% [XT,XEL,rwlset,kboth]=rwrename(XT,XEL,rwlset,editmode)
% Last modified 5-25-03
%
% Rename, or assign new ids to, one or more ring-width series; function called from rwmeas
%
%*** INPUT 
%
% XT, XEL: structures with ring-width data and associated information (see rwmeas)
%
% rwlset.  structure of .rwl file sets, with fields:
%   .name {j}s    short name or code of rwlset (e.g., padwt1}
%   .describe{j}   cell array of text description of the jth .rwl set
%   .trimall{j}    (1 x 2)i specified time coverage of data in the rwlset (all series truncated to this)
%   .idnames{j}  cell array of id names in jth set
%   .trimeach{j}  cell array of specified start and end years of spans for indiv series in the rwlset
%
% editmode (s) tells whether series to be selected from total-width set or partial-width set ('Total' or 'Partial')
% id1{} cell array of total width (XT) names currently in the storage file
% id2{} cell array of EWW (also LWW) names currently in the storage file
%
%
%*** OUTPUT
%
% Revised XT,XEL, rwlset
% kboth==1 if have answered yes to renaming the ids of the other data-type if the name exist
% kboth==0 if have answered no to that
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%


% Store series ids
id1=XT.id;
id2=XEL.id;

% Numbers of series currently stored
num1=length(id1); % number of XT ids
if num1==1 & isempty(id1{1});
    num1=0;
    if strcmp(editmode,'Total');
        error('No measured total-width series in the file; cannot delete any');
    end;
end;
num2=length(id2); % number of XEL ids
if num2==1 & isempty(id2{1});
    num2=0;
    if strcmp(editmode,'Partial');
        error('No partial-width series in the file to begin with; cannot delete any');
    end;
end;

% -- HOW MANY RWLSETS NOW EXIST

if isempty(rwlset.name{1});
    nsets=0;
else;
    nsets=length(rwlset.name);
end;

if strcmp(editmode,'Total');
    idset=id1;
    idother=id2;
    numkey=num1;
    numother=num2;
elseif strcmp(editmode,'Partial');
    idset=id2;
    idother=id1;
    numkey=num2;
    numother=num1;
else;
    error('Invalid editmode');
end;
numset=length(idset);


% % Build default menus
% yesx = cellstr(repmat('-Y',numset,1));
% nox  = cellstr(repmat('-N',numset,1));
% grpyx = cellstr([char(idset) char(yesx)]);
% grpnx = cellstr([char(idset) char(nox)]);


khow=menu('Choose','Rename series individually','Substitute 3-letter codes');


if khow==2; % substitute one 3-letter code for another
    if strcmp(editmode,'Total');
        kboth=questdlg('Apply the same substitutions to any partial-width names');
    elseif strcmp(editmode,'Partial');
        kboth=questdlg('Apply the same substitutions to any measured-total-width names');
    end;
    if strcmp(kboth,'Yes');
        kboth=1;
    else;
        kboth=0;
    end;
    stridset=char(idset);
    s1list = unique(cellstr(stridset(:,1:3))); % unique 3-letter name-starts
    nums1 = length(s1list);
    for n =1:nums1; % loop over unique 3-letter codes
        s1this=upper(s1list{n});
        kwh4=1;
        while kwh4; % While loop to allow correction of invalid substitution code
            krep = questdlg(['Substitute for the 3-letter code ' s1this '?']);
            if strcmp(krep,'Yes');
                prompt={'Enter replacement text (3 letters):'};
                def={s1this};
                dlgTitle=['Replacement code for ' s1this];
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                s1new=upper(char(answer{1}));
                kcheck=questdlg(['Substitute ' s1new ' for ' s1this '?']);
                % Build flag for new code matching some existing code
                kdupe = strmatch(s1new,s1list,'exact');
                
                if strcmp(kcheck,'Yes') & length(s1new)==3  & isletter(s1new(1)) & isempty(kdupe);
                    kwh4=0;
                    idset=strrep(idset,s1this,s1new); % replace code in ids
                    
                    % Replace code in id lists of rwlsets
                    if kboth==1 & numother>0; 
                        idother=strrep(idother,s1this,s1new);
                    else;
                    end;
                    if nsets ~=0;  % if at least one existing rwlset
                        for m = 1:nsets;
                            thesenames=rwlset.idnames{m};
                            thesenames=strrep(thesenames,s1this,s1new);
                            rwlset.idnames{m}=thesenames;
                        end;
                        clear thesenames;
                    end; %  if nsets ~=0;
                else;
                    kwh4=1;
                    if length(s1new) ~=3 |  ~isletter(s1new(1)) | ~isempty(kdupe);
                        uiwait(msgbox([s1new ' is not 3 chars beginning with a letter, or matches an existing 3-letter code'],'Try again!','modal'));
                    end;
                end; % if strcmp(kcheck,'Yes');
            else; % strrep=='No': do not replace this 3-letter code
            end; %  if strcmp(krep,'Yes');
        end; %  while kwh4;
    end; % for n =1:nums1;
    
    
else; % khow==1, rename series individually
    if strcmp(editmode,'Total');
        kboth=questdlg('Rename matching partial-width ids if you rename total-width ids?');
    elseif strcmp(editmode,'Partial');
        kboth=questdlg('Rename matching total-width ids if you rename partial-width ids?');
    end;
    if strcmp(kboth,'Yes');
        kboth=1;
    else;
        kboth=0;
    end;
    nallow=[1 1];  % allow to choose minimum of 1 series, maximum of 1 series, at a time
    yesx = cellstr(repmat('-Y',numset,1));
    nox  = cellstr(repmat('-N',numset,1));
    grpyx = cellstr([char(idset) char(yesx)]);
    grpnx = cellstr([char(idset) char(nox)]);
    strmenu='Toggle N to Y for series to be renamed';       
    C1=idset; % cell array of names
    C2=nox;
    Lin = logical(zeros(numset,1));
    [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu); % C3 is cell array of selected series
    ithis =find(Lout);
    s1this=char(C3);
    
    kwh4=1;
    while kwh4; % While loop to allow undoable substitution of an id 
        krep = questdlg(['Substitute a new id name for  ' char(C3) '?']);
        if strcmp(krep,'Yes');
            prompt={'Enter replacement idname:'};
            def={s1this};
            dlgTitle=['New id to replace ' s1this];
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            s1new=upper(char(answer{1}));
            kcheck=questdlg(['Rename ' s1this ' as ' s1new '?']);
            
            kdupe=strmatch(s1new,idset,'exact');
            if strcmp(kcheck,'Yes') & length(s1new)>=9;
                 kwh4=1;
                 uiwait(msgbox([s1new ' is length 9 or greater'],'Try again!','modal'));
             elseif strcmp(kcheck,'Yes') & ~isletter(s1new(1));
                 kwh4=1;
                 uiwait(msgbox([s1new ' does not begin with a letter'],'Try again!','modal'));
             elseif ~isempty(kdupe);
                kwh4=1;
                uiwait(msgbox([s1new ' is already used as an idname '],'Try again!','modal'));
            else;
                kwh4=0;
                idset=strrep(idset,s1this,s1new); % rename id
                
                % Optionally rename "sister" series in other data type 
                if kboth==1 & numother>0; 
                    idother=strrep(idother,s1this,s1new);
                    if ~isempty(strmatch(s1new,idother,'exact'));
                        error(['The other type of data already has an id ' s1new]);
                    end;
                else;
                end;
                
                % HANDLE any RWSETS
                if nsets ~=0;  % if at least one existing rwlset
                    for m = 1:nsets;
                        thesenames=rwlset.idnames{m};
                        thesenames=strrep(thesenames,s1this,s1new);
                        rwlset.idnames{m}=thesenames;
                    end;
                    clear thesenames;
                end; %  if nsets ~=0;
                
            end; % if strcmp(kcheck,'Yes');
            
            
            
            
        else; % strrep=='No': do not replace this 3-letter code
        end; %  if strcmp(krep,'Yes');
    end; %  while kwh4;

end; % if khow==2; % substitute one 3-letter code for another

% Put the revised ids in XT.id 
if strcmp(editmode,'Total');
    XT.id=idset;
    if kboth==1 & numother>0;
        XEL.id=idother;
    end;
elseif strcmp(editmode,'Partial');
    XEL.id=idset;
    if kboth==1 & numother>0;
        XT.id=idother;
    end;
end;


% If no "sister" data, no need to have kboth flag set
if numother==0;
    kboth=0;
end;






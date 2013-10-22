function [X,rwlset]=rwdelete(X,rwlset,editmode)
% rwdelete:  delete measureement series from an rwmeas measurement file
% [X,rwlset]=rwdelete(X,rwlset,editmode);
% Last modified 5-26-03
%
% Delete measurement series from an rwmeas measurement file
%
%*** INPUT 
%
% X: structure with ring-width data and associated information (see rwmeas)
%
% rwlset.  structure of .rwl file sets, with fields:
%   .name {j}s    short name or code of rwlset (e.g., padwt1}
%   .describe{j}   cell array of text description of the jth .rwl set
%   .trimall{j}    (1 x 2)i specified time coverage of data in the rwlset (all series truncated to this)
%   .idnames{j}  cell array of id names in jth set
%   .trimeach{j}  cell array of specified start and end years of spans for indiv series in the rwlset
%
% editmode (s) tells whether series to be selected from total-width set or partial-width set ('Total' or 'Partial')
%
%
%*** OUTPUT
%
% Revised X, rwlset
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% rwlset not updated to reflect changed makeup of ids within this function


% Store series ids
id1=X.id;
idset=id1; % moot renaming

% Numbers of series currently stored
num1=length(id1); % number of  ids
if num1==1 & isempty(id1{1});
    num1=0;
    if strcmp(editmode,'Total');
        error('No measured total-width series in the file; cannot delete any');
    elseif strcmp(editmode,'Partial');
         error('No measured total-width series in the file; cannot delete any');
    end;
end;

% Store data
D=X.data;

% -- HOW MANY RWLSETS NOW EXIST
if isempty(rwlset.name{1});
    nsets=0;
else;
    nsets=length(rwlset.name);
end;

numids=length(idset); % number of series (ids)
khow=menu('Choose',['Select individual ' lower(editmode) '-width series to delete'],['Delete all ' lower(editmode) '-width series']);

if khow==2; % Delete all series
    
    kwh1=1;
    while kwh1;
        kreally=questdlg(['Really want to delete all ' lower(editmode) '-width series?']);
        if strcmp(kreally,'Yes');
            X.id={[]};
            X.data={[]};
            X.span={[]};
            X.who={[]};
            X.when={[]};
            X.remeasure={[]};
            X.summary=[];
            kwh1=0;
            strmess1={['Although you are deleting all ' lower(editmode) '-width series,'],...
                        'any existing rwlsets will not be changed.  Thus after this operation, ',...
                        'an rwlset could conceivably contain the id of a series deleted from the',...
                        'measurement file.  You would need to edit the rwlset to fix that'};
            uiwait(msgbox(strmess1,'Message','modal'));
        else;
            clc;
            return;
        end; % if strcmp(kreally,'Yes');
    end; % while kwh1;
    
else; % khow==1, delete series individually
    
    nallow=[0 (numids-1)];  % allow to choose minimum of 1 series, maximum of all but 1 series
    yesx = cellstr(repmat('-Y',numids,1));
    nox  = cellstr(repmat('-N',numids,1));
    grpyx = cellstr([char(idset) char(yesx)]);
    grpnx = cellstr([char(idset) char(nox)]);
    strmenu='Toggle N to Y for series to be deleted';       
    C1=idset; % cell array of names
    C2=nox;
    Lin = logical(zeros(numids,1));
    kwh2=1;
    while kwh2;
        [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu); % C3 is cell array of selected series
        iout =find(Lout);
        
        if isempty(iout); % no series selected
            kquest1=questdlg('Do you really mean to delete no series?');
            if strcmp(kquest1,'Yes');
                kwh2=0;
                % X and rwlset remain as is
                % no action needed
            else;
                kwh2=1
            end;
        else; % some series were selected
            strlist = C3;
            uiwait(msgbox(strlist,'To Delete:','modal'));
            kkill=questdlg('Proceed with deleting those series?');
            if strcmp(kkill,'Yes');  
                kwh2=0;
                X.id(iout)=[];
                X.data(iout)=[];
                X.span(iout)=[];
                X.who(iout)=[];
                X.when(iout)=[];
                X.remeasure(iout)=[];
                              
            else; %  if strcmp(kkill,'Yes');
                kwh2=1;
            end; % if strcmp(kkill,'Yes');
        end; 
        
    end; % while kwh2;
end; % if khow==2; 

    
   
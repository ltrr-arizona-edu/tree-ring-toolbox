function [id1,d1,id2,d2]=cullids(id1,d1,id2,d2,rwlset,ftype)
% cullids:  select subset of ring-width ids using an rwlset
% [id1,d1,id2,d2]=cullids(id1,d1,id2,d2,rwlset,ftype);
% Last revised 2-11-02
%
% Select subset of ring-width ids using an rwlset.  Function called by rwc2rwl.m that
% allows user to pick apply a filter to the ids in a ring-width set. The filter is a 
% previously built rwlset.  This utility function part of the rwmeas suite.
%
%*** INPUT - see notes
% 
% id1{} cell array of ids for total-width series
% d1{} cell array of 2-col time series matrices  (year as col 1) for series in id1
% id2{} cell array of ids for total-width series
% d2{} cell array of 2-col time series matrices  (year as col 1) for series in id2
% rwlset{} structure with information on data-set specifications for rwl groups of series
% ftype (1 x ?)s  type of data to cull ids for; may be 'ew', 'lw', 'rwc' or 'rwm'
%
%
%*** OUTPUT
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% See calling function rwc2rwl.m for information on id1, d1, id2, d2
% See rwmeas.m for information on rwlset and ftype

%--- check inputs

% ftype
if ~ischar(ftype);
    error('ftype must be string');
end;
if ~any(strcmp({'ew','lw','rw','rwc','rwm'},ftype));
    error('ftype must be ew, lw rwc or rwm');
end;

% rwlset
if ~isstruct(rwlset);
    error('rwlset not a structure');
end;
if length(rwlset.name)==1 & isempty(rwlset.name{1});
    error('rwlset is empty');
else;
    nsets = length(rwlset.name);
end;

%--- CHOOSE THE RWLSET TO USE

C1=(rwlset.name)';
C2=cellstr(repmat('-N',nsets,1));
Lin = logical(zeros(nsets,1));
nallow=[1 1];
strmenu='Toggle the desired rwlset to Y';
[C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu);
% C3 is now the name of the set;  Lout is pointer to the set


%--- Get desired ids for the set
idwant = rwlset.idnames{Lout}; % cell of strings
nwant = length(idwant); % how many desired series in the set

%--- Get the desired periods for series in the set
pereach = rwlset.trimeach{Lout}; % cell of periods

%-- Get the overall restrictions
perall = rwlset.trimall{Lout};

%--- Compute the start and end year of restricted output. Trimall trumps trimeach
T = repmat(NaN,nwant,2);
for n =1:nwant;
    per1 = pereach{n};
    if isempty (per1) & isempty(perall);
        T(n,:)=[-inf inf];
    elseif isempty(per1); % no individual restrictions, but a global
        T(n,:)=perall;
    elseif isempty(perall); % no global, but an individual
        T(n,:)=per1;
    else; % both global and indiv
        T(n,:)=[max([per1(1) perall(1)]) min([per1(2) perall(2)])];
    end;
end;
% Message if truncation of period eliminates any series
delta1 = (diff(T'))';
if any(delta1<0);
    jcut=find(delta1<0);
    mess1='Trimmed period eliminates one or more series from the rwlset writeout';
    uiwait(msgbox(mess1,'Message','modal'));
else;
    jcut=[];
end;


%--- Trim id1, id2, d1, d2 to contain only elements for the set

nid1=length(id1);
if nid1==1 & isempty(id1{1});
    n1=0;
    % no change to id1 or d1
else;
    Lkeep = logical(zeros(nid1,1));
    for n=1:nwant;
        Tthis=T(n,:);
        idthis=idwant{n};
        L = strcmp(id1,idthis);
        iseq=find(L); % seq number of idthis in id1
        Lkeep(L)=1;
        if ~isempty(iseq);
            % Trim data in d1
            x=d1{iseq};
            Lgood = x(:,1)>=Tthis(1) & x(:,1)<=Tthis(2);
            x=x(Lgood,:);
            d1{iseq}=x;
            if ~isempty(jcut);
                if any(n==jcut);
                    Lkeep(L)=0;
                end;
            end;
        end;
    end; 
    id1=id1(Lkeep);
    d1=d1(Lkeep);
    n1=sum(Lkeep);
end;

nid2=length(id2);
if nid2==1 & isempty(id2{1});
    n2=0;
    % no change to id2 or d2
else;
    Lkeep = logical(zeros(nid2,1));
    for n=1:nwant;
        Tthis=T(n,:);
        idthis=idwant{n};
        L = strcmp(id2,idthis);
        iseq=find(L); % % seq number of idthis in id2
        Lkeep(L)=1;
        % Trim data in d2, if the series is in id2
        if ~isempty(iseq);
            x=d2{iseq};
            Lgood = x(:,1)>=Tthis(1) & x(:,1)<=Tthis(2);
            x=x(Lgood,:);
            d2{iseq}=x;
            
            if ~isempty(jcut);
                if any(n==jcut);
                    Lkeep(L)=0;
                end;
            end;
        end;
    end; 
    id2=id2(Lkeep);
    d2=d2(Lkeep);
    n2=sum(Lkeep);
end;


if any(strcmp({'ew','lw','rwc'},ftype))  & n2==0;
     error('No overlap of rwlset with data of ew/lw');
 end;
if any(strcmp({'rwm'},ftype))  & n1==0;
     error('No overlap of rwlset with data of total width');
end;
if n1==0 & n2==0;
    error('No series match ids in rwlset');
end;



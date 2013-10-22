function X=rwimport(X,editmode,kdigits)
%
% Last Revised 2006-08-09
%NEXT
% * code for culling ids depending on what already in X.id and whether replace only, add only, or add and replace mode
% * above that, code for select-series inport
% see *** IF SELECTED

% rwimport: import .rw, .eww or .lww data into cell-format ringwidth storage
% X=rwimport(X,editmode);
% Last revised 10-8-01
%
% Import .rw, .eww or .lww data into cell-format ring-width storage.  You want to incorporate measurements
% made elsewhere into a cell-format file.
%
%*** INPUT 
%
% X -- structure of total or partial ring width data (see rwmeas.m), or [] if no such structure yet
%   exists
% editmode (1 x ?)s type of measurement being imported
%       =='Total'
%       =='Partial'
% kdigits   ==1 input data units 100ths of mm
%           ==2 input data units 1000ths of mm
%
%
%*** OUTPUT
% 
% X -- revised structure
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED 
%
% dirfls
%
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Rev2006-08-09: new in arg kdigits. Used in call to read5 to allow 1000ths mm as well as 100ths mm input units

%--- COMPUTE CURRENT NUMBER OF STORED SERIES

if isempty(X.id{1});
    nser=0;
else;
    nser = length(X.id);
end;
mnext = nser+1;

%***  SET NUMBER OF FILE-TYPES IMPORT

switch editmode;
case 'Total';
    ntype=1; % just .rw
    suff={'rw'};
case 'Partial';
    ntype=2; % .eww and .lww
    suff={'eww','lww'};
otherwise;
end;


%*** UI TO GET PATH OF IMPORT FILES (DEFAULT C:\RWIMPORT\ -->path1

kmen1=menu('Choose one',...
    'Import files from c:\import\',...
    'Browse and click to get path for import');
if kmen1==1;
    path1='c:\import\';

else;
    if strcmp(editmode,'Total');
        [file1,path1]=uigetfile('*.rw','Click on any file in desired input directory');
    elseif strcmp(editmode,'Partial');
        [file1,path1]=uigetfile({'*.eww';'*.lww'},'Click on any file in desired input directory');
    end;
    clear file1;
end;



%*** GET THE PREFIXES OF FILENAMES , AND CHECK THAT AT LEAST ONE OF REQUIRED FILE TYPE EIISTS IN DIRECTORY; SORT

pre=dirfls4(path1,suff);
gre = pre; % store the original versions of names 
%  If total width, pre.x1{} now holds filenames
% If partial, pre.x1{} holds .eww names and pre.x2{} holds latewood names

% Check that some files of each specified type have been found
for n=1:ntype;
    eval(['xtemp=pre.x' int2str(n) ';']);
    if isempty(xtemp);
        error(['No files with suffix ' suff{1} ' in ' path1]);
    else;
    end;
end;
clear xtemp;

% If partial, check that number of .eww files is same as number of .lww files in the import directory
if ntype==2;
    if length(pre.x1) ~= length(pre.x2);
        error([num2str(length(pre.x1)) ' .eww files, but  ' num2str(length(pre.x2)) ' .lww files in ' path1]);
    end;
end;

% Build cross-ref matrix of ids to filenames and store sorted ids in pre.x1 and pre.x2
for n=1:ntype;
    eval(['nms=char(pre.x' int2str(n) ');']);
    nms_orig=cellstr(nms); 
    nms1=fixname1(nms);
    eval(['pre.x' int2str(n) '=cellstr(nms1);']);
    eval(['gre.x' int2str(n) '=cellstr(nms);']); % original names, before fixnam
    eval(['FX' int2str(n) '= [pre.x' int2str(n)  ' cellstr(nms) ];']);
    eval(['xthis=pre.x' int2str(n) ';']);
    
    % Sort the "corrected" names, but must also keep matching list of uncorrected names because that is what
    % the filew (rw, etc) are stored as
    [xthis,i1]=sort(xthis);  % sort the series names
    eval(['pre.x' int2str(n) '=xthis;']);
    ythis = nms_orig(i1);
    eval(['gre.x' int2str(n) '=ythis;']);
end;
clear xthis nms nms1 n ;


% IF PARTIAL, CHECK THAT FILES OF IDENTICAL PREFIX FOR EWW AND LWWW
nfiles = length(pre.x1);
if strcmp(editmode,'Partial');
    if ~all(strcmpi(pre.x1,pre.x2));
        disp(char(pre.x1));
        disp(char(pre.x2));
        error(['Look!  .eww names not identical to .lww names']);
    end;
end;



%*** CHECK THAT FILENAMES VALID CODING OF SITE/TREE/CORE/SEQ;  

% Run corenm1.  Will bomb with error message if any id invalid
id = pre.x1; % will send either the .rw or .eww names (.lww are identical to .eww)
corenm1(id);


%*** QUERY FOR BULK IMPORT VS SELECTED-FILE IMPORT

kmen=menu('Choose import mode',...
    'Bulk -- all files of selected type in the input directory',...
    'Selected -- will menu-pick those files to import');
if kmen==1;
    howin ='Bulk';
else;
    howin='Select';
end;
clear kmen;


%*** QUERY FOR OVERWRITE VS ADD ONLY (OR SOME COMBINATION)

kmen=menu('Choose one',...
    'Add or overwrite',...
    'Add only',...
    'Overwrite only');
if kmen==1;
    ksub = 'AddOver';
elseif kmen==2;
    ksub='AddOnly';
elseif kmen==3;
    ksub='OverOnly';
end;


%*** IF SELECTED IMPORT, CULL id.x1

if strcmp(howin,'Bulk');
elseif strcmp(howin,'Select');
    Lok = logical(zeros(nfiles,1));
    flagno=repmat('-N',length(pre.x1),1); % string of Nos for selection
    pickme = char([char(pre.x1)  flagno],'Accept and Quit');
    pickme=cellstr(pickme);
    kwh1=1;
    while kwh1==1;
        kmen = menu('Toggle to ''y'' to include',pickme);
        if kmen==nfiles+1;
            if all(Lok==0);
                uiwait(msgbox('You Still Have Not Picked Any Series','Message','modal'));
            else;
                kwh1=0;
            end;
        else;
            thisone =flagno(kmen,:);
            if strcmp(thisone,'-N');
                thisone='-Y';
                Lok(kmen)=1;
            elseif strcmp(thisone,'-Y');
                thisone='-N';
                Lok(kmen)=0;
            end;
            flagno(kmen,:)=thisone;
            pickme = char([char(pre.x1)  flagno],'Accept and Quit');
            pickme=cellstr(pickme);
        end;
    end;
    
    pre.x1=pre.x1(Lok);
    if ntype==2;
        pre.x2=pre.x2(Lok);
    end;
    nfiles=length(pre.x1); % recompute number of files to consider importing
end;




%*** CULL pre.x1, pre.x2, nfiles depending on makeup of X.id and mode for replace/add

Lkill = logical(zeros(nfiles,1)); % pre.x1 and pre.x2 elements to delete from consideration
if strcmp(ksub,'OverOnly');    % Replace only: check pre.x1 vs X.id; cull only pre.x1 in X.id
    if isempty(X.id{1}); % No series yet stored
        error('Replace-only mode, but no series initially there to replace');
    else; % X.id{1} not empty
        for n=1:nfiles; % loop over selected files
            idthis = pre.x1{n};
            if ~any(strcmpi(idthis,X.id));
                Lkill(n)=1; % flag for delete
            end;
        end;
    end;
elseif strcmp(ksub,'AddOnly'); % Add only .. ditto, but inverse
    if isempty(X.id{1}); % No series yet stored.  So will accept all of pre.x1
        % No action needed;  Lkill remains all zero
    else; % X.id{1} not empty
        for n=1:nfiles; % loop over selected files
            idthis = pre.x1{n};
            if any(strcmpi(idthis,X.id));
                Lkill(n)=1; % flag for delete
            end;
        end;
    end;
elseif strcmp(ksub,'AddOver'); % Add or replace -- use all id.x1
    % No action needed;  accept Lkill as all-zero
else;
    error('Invalid ksub string');
end;

% Cull the selected ids

if ~any(Lkill);  % no reduction of pre.x1 and pre.x2 needed
else;
    pre.x1(Lkill)=[];
    pre.x2(Lkill)=[];
    nfiles = length(pre.x1);
end;
if nfiles==0;
    error('You have ended up with NO series to be added or replaced');
end;
    



%*** IMPORT THE DATA INTO X

if strcmp(editmode,'Total');
    for n =1:nfiles;
        id=pre.x1{n};
        id_orig =gre.x1{n};
        
        % *******Check id against current ids & set index for storage
        
        m=find(strcmpi(id,X.id));
        if isempty(m);
            m=mnext;
            mnext=mnext+1;
        else;
            m=m;
        end;
                
        pfinput = [path1 id_orig '.' suff{1}];
        [x,person,when]=rwread5(pfinput,kdigits);
        
        
        % Store 
       
        X.data{m}=x;
        X.who{m}=person;
        X.when{m}=when;
        X.id{m}=id;
        X.remeasure{m}=[];
        
        disp([pfinput ' imported']);
    end;
else; % editmode Partial
    
    for n =1:nfiles;
        id=pre.x1{n};
        id_orig=gre.x1{n};
        
        % *******Check id against current ids & set index for storage
        
        m=find(strcmpi(id,X.id));
        if isempty(m);
            m=mnext;
            mnext=mnext+1;
        else;
            m=m;
        end;
        
        
        % Earlywood
        iwinner=strmatch(id,FX1(:,1)); % find row of this series in file cross-reference matrix
        idwinner = FX1(iwinner,2); % get the original filename
        pfinput = [path1 char(idwinner) '.' suff{1}];
        [x1,person,when]=rwread5(pfinput,kdigits);
        
        % Latewood
        iwinner=strmatch(id,FX2(:,1)); % find row of this series in file cross-reference matrix
        idwinner = FX2(iwinner,2); % get the original filename
        pfinput = [path1 char(idwinner) '.' suff{2}];
        [x2,person,when]=rwread5(pfinput,kdigits);
        
        % Check start year
        if ~ (  x1(1,1)==x2(1,1) | x1(1,1)==(x2(1,1)+1)  );  
            strwarn=[pfinput ' start yrs of eww and lww : ' int2str([x1(1,1) x2(1,1)])];
            uiwait(msgbox(strwarn,'Fatal Error','modal'));
            error([id ': inconsistent first years of early and latewood']);
        end;
        % Check end year
        if ~ (  x1(end,1)==x2(end,1) | x1(end,1)==(x2(end,1)+1)  );  
            strwarn=[pfinput ' end yrs of eww and lww : ' int2str([x1(end,1) x2(end,1)])];
            uiwait(msgbox(strwarn,'Fatal Error','modal'));
            error([id ': inconsistent last years of early and latewood']);
        end;
        
        % Compute start year for tsm
        yrgo = min([x1(1,1)  x2(1,1)]);
        yrsp = max([x1(end,1)  x2(end,1)]);
        yrx = (yrgo:yrsp)';
        nyrx = yrsp-yrgo+1;
        
        % Reserve stoage
        x=repmat(NaN,nyrx,4);
        x(:,1)=yrx;
        
        % store early
        nsize=size(x1,1);
        i1=x1(1)-yrgo+1;
        i2=i1+nsize-1;
        x(i1:i2,2)=x1(:,2);
        
        % store late
        nsize=size(x2,1);
        i1=x2(1)-yrgo+1;
        i2=i1+nsize-1;
        x(i1:i2,3)=x2(:,2);
        
        % computed total
        x(:,4)=x(:,2)+x(:,3);
        
                
        % Store 
       
        X.data{m}=x;
        X.who{m}=person;
        X.when{m}=when;
        X.id{m}=id;
        X.remeasure{m}=[];
        
    end;
end;

% Inform how many series imported (overwritten and added)
strmess1=[int2str(nfiles)  ' .rw files successfully imported by rwimport.m'];
strmess2=[int2str(nfiles)  ' .eww and ' int2str(nfiles) ' .lww files  successfully imported by rwimport.m'];
if strcmp(editmode,'Total');
    strmess=strmess1;
else;
    strmess=strmess2;
end;
uiwait(msgbox(strmess,'Message','modal'));


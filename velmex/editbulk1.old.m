function X=editbulk(X);
% editbulk: bulk edit some feature of ringwidth
% X=editbulk(X,bulkmode);
% Last Revised 10-01-01
%
% Apply same change to multiple series: c hange the last year, insert a locally absent ring, 
% or remove a locally absent ring
%
%
%*** INPUT
%
% X -- structure of total or partial width (see rwmeas)
%
%*** OUTPUT
%
% X -- revised structure
%
%*** NOTES
%
% Use can select to 1) edit all series, 2) select those to edit from click window, 3) select those not to be edited from 
% click window
%
% No series include for bulk edit us allowed to have any remeasurement pending (i.e., the X.remeasure{} 
% for series to be bulk-edited must be [].

%--- HAVE DATA?

id =X.id;
d=X.data; % data
p=X.who; % person
w=X.when; % date
dspan =X.span; % span
remeas=X.remeasure;
if isempty(d{1});
    uiwait(msgbox('X EMPTY -- CANNOT EDIT','Message','modal'));
    return;
end;


%--- WHAT CHANGE

kchng=menu('Choose ',...
    'CHANGE END YEAR',...
    'INSERT LOCALLY ABSENT RING, KEEPING SAME START YEAR',...
    'INSERT LOCALLY ABSENT RING, KEEPING SAME END YEAR',...
    'REMOVE LOCALLY ABSENT RING, KEEPING SAME START YEAR',...
    'REMOVE LOCALLY ABSENT RING, KEEPING SAME END YEAR');
if kchng==1; % CHANGE END YEAR
   prompt={'Enter desired end year:'};
   def={'2001'};
   dlgTitle='Set the End Year for All Selected Series';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
   yrend = str2num(answer{1});
   yrnow = datevec(date);
   yrnow=yrnow(1);
   if yrend>yrnow;
       error([num2str(yrend) ' is a future year']);
   end;
elseif kchng==2 | kchng==3;; % INSERT LOCALLY ABSENT RING
   prompt={'Enter year:'};
   def={'2001'};
   dlgTitle='Insert LA ring AFTER  specified year';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
   yrkey = str2num(answer{1});
   yrnow = datevec(date);
   yrnow=yrnow(1);
   if yrkey>=yrnow;
       error([num2str(yrkey) ' is current or future year']);
   end;
elseif kchng==4 | kchng==5;; % REMOVE LOCALLY ABSENT RING
   prompt={'Enter year:'};
   def={'2001'};
   dlgTitle='Remove LA ring currently in specified year?';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
   yrkey = str2num(answer{1});
   yrnow = datevec(date);
   yrnow=yrnow(1);
   if yrkey>=yrnow;
       error([num2str(yrkey) ' is current or future year']);
   end; 
end;


%----- EDIT WHAT TYPE OF MEASUREMENT?

if size(d{1},2)==2;
    editmode='Total';
    lab1='RING WIDTH (mm)';
    jcol=1; % data col 
    jj=1;
elseif size(d{1},2)==4;
    editmode='Partial';
    lab1='EWW+LWW (mm)';
    jcol =2; % data col, total width
    jj=[1 2 3];
end;


%---- GET INDICES TO SERIES TO BE EDITED

% Build menu choices
nser1 = length(id); % total number of available series
idlist=id;
idlist{nser1+1}='FINISHED';
yesman = cellstr(repmat('-Y',nser1+1,1));
noman = cellstr(repmat('-N',nser1+1,1));
yesman{nser1+1}='';
noman{nser1+1}='';
idyes = cellstr([char(idlist) char(yesman)]);
idno = cellstr([char(idlist) char(noman)]);

kmen1=menu('Choose one',...
    'APPLY CHANGE TO ALL SERIES',...
    'CLICK ON SERIES TO BE CHANGED',...
    'CLICK ON SERIES TO BE EXCLUDED FROM CHANGE');
if kmen1==1; % SAME CHANGE TO ALL SERIES
    jedit=(1:nser1)'; % edit these
elseif kmen1==2; % CLICK ON SERIES TO CHANGE
    C=idno;
    N=noman;
    kwh=1; % while control
    L=logical(zeros(nser1+1,1));
    while kwh==1;
        kmen2 = menu('Click so that series to be changed marked with -Y',C);
        if kmen2==nser1+1;
            jedit=find(L);
            kwh=0;
        else; 
            if strcmp(N{kmen2},'-N');
                N{kmen2}='-Y';
                L(kmen2)=1;
                if ~isempty(remeas{kmen2});
                    uiwait(msgbox([id{kmen2} 'flagged for remeasuring -- cannot bulk edit it'],'Invalid Choice','modal'));
                    L(kmen2)=0;
                    N{kmen2}='-N';
                end;
            elseif strcmp(N{kmen2},'-Y');
                N{kmen2}='-N';
                L(kmen2)=0;
            end;
            C = cellstr([char(idlist) char(N)]);
        end;
    end;
elseif kmen1==3; % CLICK ON SERIES TO EXCLUDE FROM CHANGE
    C=idyes;
    N=yesman;
    kwh=1; % while control
    L=logical(ones(nser1+1,1));
    while kwh==1;
        kmen2 = menu('Click so that series to be changed marked with -Y',C);
        if kmen2==nser1+1;
            jedit=find(L);
            kwh=0;
        else;
            if strcmp(N{kmen2},'-N');
                N{kmen2}='-Y';
                L(kmen2)=1;
                if ~isempty(remeas{kmen2});
                    uiwait(msgbox([id{kmen2} 'flagged for remeasuring -- cannot bulk edit it'],'Invalid Choice','modal'));
                    L(kmen2)=0;
                    N{kmen2}='-N';
                end;
            elseif strcmp(N{kmen2},'-Y');
                N{kmen2}='-N';
                L(kmen2)=0;
            end;
            C = cellstr([char(idlist) char(N)]);
        end;
    end;
    
end; % kmen1=menu('Choose one',...


%--- CHECK TO BE SURE NO SELECTED SERIES ARE FLAGGED FOR REMEASURING

if isempty(jedit);
    uiwait(msgbox('No series have been marked for change','Message','modal')); 
    return;
end;


%--  CHECK TO SEE THAT SELECTED YEAR IS COMPATIBLE WITH YEAR COVERAGE AND DATA OF SERIES

npick = length(jedit);
yrs = repmat(NaN,npick,2);
if kchng==2 | kchng==3; % Insert LA ring
    for n=1:npick;
        jthis=jedit(n);
        yrs (n,1)=X.span{jthis}(1);
        yrs (n,2)=X.span{jthis}(2);
    end;
    if any(yrkey<yrs(:,1))  | any(yrkey>=yrs(:,2));
        error(['Year ' num2str(yrkey) ' out of range for at least one selected series']);
    end;
elseif kchng==4 | kchng==5; % remove LA ring
    xkey = repmat(NaN,npick,1);  % to store total ring-width for selected year 
    for n=1:npick;
        jthis=jedit(n);
        yrs (n,1)=X.span{jthis}(1);
        yrs (n,2)=X.span{jthis}(2);
        ikey = d{jthis}(:,1)==yrkey;
        if  ~any(ikey);
            error([id{jthis} ' does not have data for ' num2str(yrkey)]);
        end;
        if strcmp(editmode,'Total');
            xkey(n) = d{jthis}(ikey,2);
        elseif strcmp(editmode,'Partial');
            xkey(n) = d{jthis}(ikey,4);
        end;
    end;
    if any(yrkey<yrs(:,1))  | any(yrkey>=yrs(:,2));
        error([num2str(yrkey) ' outside range of some selected series']);
    end;
    if ~all(xkey==0);
        error(['Not all selected series have a LA ring for year ' num2str(yrkey)]);
    end;
end;

%-- MAKE THE CHANGE

for n = 1:npick;
    jthis=jedit(n);
    
    if strcmp(editmode,'Total');
        jcol = 2;
        xzero=0;
    elseif strcmp(editmode,'Partial');
        jcol = [2 3 4];
        xzero=[0 0 0];
    end;
    dthis = d{jthis}; % data -- 2 cols or 4 cols
    yrthis= dthis(:,1);
    yron1 = yrthis(1);
    yroff1= yrthis(end);
    nthis =length(yrthis);
    i1 =   flipud((0:(nthis-1))');
    sthis = dspan{jthis}; % span
    
    if kchng==1; % change last year
        yrnew = yrend-i1; % revised year vector;
        dthis(:,1)=yrnew;
        d{jthis}=dthis;
    elseif kchng==2 | kchng==3; % INSERT LA RING,
        i1 = flipud((0:(nthis))');
        ihere = yrthis==yrkey;
        ibefore = yrthis<=yrkey;
        iafter = yrthis>yrkey;
        xblock1 = dthis(ibefore,jcol);
        xblock2 = [ xzero];
        xblock3= dthis(iafter,jcol);
        xblock = [xblock1; xblock2; xblock3];
        if kchng==2; % KEEP START YEAR
            yrend=yroff1+1;
        elseif kchng==3; % KEEP END YEAR 
            yrend=yroff1;
        end;
        yrnew = yrend-i1;
        d{jthis}=[yrnew xblock];
    elseif kchng==4 | kchng==5; % REMOVE LA RING
        i1 = flipud((0:(nthis-2))');
        ihere = yrthis==yrkey;
        ibefore = yrthis<yrkey;
        iafter = yrthis>yrkey;
        xblock1 = dthis(ibefore,jcol);
        xblock3= dthis(iafter,jcol);
        xblock = [xblock1;  xblock3];
        if kchng==4; % KEEP START YEAR
            yrend=yroff1-1;
        elseif kchng==5; % KEEP END YEAR 
            yrend=yroff1;
        end;
        yrnew = yrend-i1;
        d{jthis}=[yrnew xblock];
    end;
    son = yrnew(1);
    soff = yrnew(end);
    dspan{jthis}=[son soff];
    X.span=dspan;
    X.data=d;
end;


    
    
    
    
        
        
    
    



    
        


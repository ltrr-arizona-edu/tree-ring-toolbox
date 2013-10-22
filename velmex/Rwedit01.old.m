function V=rwedit01(V)
% rwedit01:  edit cell-formatted ring-width measurements
% V=rwedit(V);
% Last revised 9-28-01
% 
% Edit cell-formatted ring-width data.  Changes allowed include:
%-change start or end year
%-add or remove a locally absent ring (0 value)
%-change a mis-identified real ring to a false ring
%-change a mis-identified false ring to a real ring
%-add or remove a contiguous stretch of missing values (NaN)
%-change a measured value
%
%
%*** INPUT
%
% V -- structure of total or partial ring-width data (see notes)
%
%*** OUTPUT 
%
% V - revised structure
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Template was rwedit from \rwlook\
% Separate partial-ring components are not editable within edit mode -- must measure
% The fields T.remeasure and P.remeasure are changed as necessary to record which years for which series should be 
%   remeasured later.  For example, if you split a real ring into two real rings (false to real), you will eventually 
%   need to go back and remeasure the total and partial portions of those affected rings.  This fuction just assigns
%   arbitrary splits.
% One change only per series is allowed in a single run of rwedit01.  If the .remeasure field is not empty, you are not 
%   allowed to edit the series.

%--- HAVE DATA?

id = V.id;
d=V.data; % data
p=V.who; % person
w=V.when; % date
dspan = V.span; % span
remeas=V.remeasure;
if isempty(d{1});
    uiwait(msgbox('V EMPTY -- CANNOT EDIT','Message','modal'));
    return;
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

%------ WHICH SERIES

kser=menu('CHOOSE SERIES TO EDIT',id);


%--- CHECK THAT ANY REQUIRED MEASUREMENTS HAVE BEEN MADE FROM PREVIOUS EDITING

if ~isempty(V.remeasure{kser});
    uiwait(msgbox([id{kser} ':  CANNOT EDIT TILL YOU REMEASURE FLAGGED YEARS'],'Message','modal'));
    return;
end;

%--- GET THE SERIES 

X0=d{kser}; % either 2-col or 4-col
yrx0 = X0(:,1); % store year
X0(:,1)=[]; % strip off year column
x0=X0(:,jcol); % store total-width column separately

X=X0;
x=x0;
yrx=yrx0;
% Plot the series
subfcn01(x0,yrx0,id{kser},lab1);

kquest=questdlg('IS THIS THE SERIES YOU WANT TO EDIT?');
switch kquest;
case 'Yes';
case {'No','Cancell'};
    return;
end;


%-- HOW TO EDIT

kwh1=1;
kdone=0;
while kwh1==1;
    kmen1 = menu('What do you want to do?',...
        'CHANGE END YEAR (SHIFT SERIES)',...
        'INSERT OR REMOVE LOCALLY ABSENT RING',...
        'CHANGE REAL RING TO FALSE',...
        'CHANGE FALSE RING TO REAL (SPLIT THE RING)',...
        'INSERT OR REMOVE A STRECH OF MISSING VALUES',...
        'CHANGE A MEASURED VALUE',...
        'ABORT THE CHANGE AND START OVER',...
        'ACCEPT THE CHANGE AND QUIT');
    
    if kmen1==1; % CHANGE END YEAR (SHIFT SERIES)
        if kdone==1;
            strmsg1='CANNOT STACK CHANGE ON CHANGE: ABORT PREVIOUS CHANGE FIRST';
            uiwait(msgbox(strmsg1,'Message','modal'));
            
        else;
            yrsp = yrx(end); 
            yrgo = yrx(1);
            syrgo = num2str(yrgo);
            syrsp = num2str(yrsp);
            txttit1=['Change start year from ' syrgo];
            txttit2=['Change end year from ' syrsp];
            kmen2=menu('Choose One',...
                'SPECIFY NEW START YEAR',...
                'SPECIFY NEW END YEAR');
            if kmen2==1;
                prompt={'Enter new start year:'};
                def={syrgo};
                dlgTitle=txttit1;
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                yron = str2num(answer{1});
                % Change data
                yrshift = yron-yrx(1); % shift in years
            elseif kmen2==2; 
                prompt={'Enter new end year:'};
                def={syrsp};
                dlgTitle=txttit2;
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                yroff = str2num(answer{1});
                % Change data
                yrshift = yroff-yrx(end); % shift in years
            end;
            clear kmen2;
            % Revised data
            yrx = yrx+yrshift;
            kdone=1;
            % Plot series
            subfcn02(x0,yrx0,x,yrx,id{kser},lab1);
            
        end; % if kdone==1;
    elseif kmen1==2; % ADD OR REMOVE LOCALLY ABSENT RING
        if kdone==1;
            strmsg1='CANNOT STACK CHANGE ON CHANGE: ABORT PREVIOUS CHANGE FIRST';
            uiwait(msgbox(strmsg1,'Message','modal'));
        else;
            yrsp = yrx(end); 
            yrgo = yrx(1);
            syrkey='9999';
            txttit1=['Locally Absent Ring Operation'];
            kmen2=menu('Locally Absent Ring',...
                'INSERT, KEEPING SAME START YEAR',...
                'INSERT, KEEPING SAME END YEAR',...
                'REMOVE, KEEPING SAME START YEAR',...
                'REMOVE, KEEPING SAME END YEAR');
            % Which year
            if kmen2==1 | kmen2==2;
                str1='INSERT';
            else;
                str1='REMOVE';
            end;
            
            % PROMPT FOR YEAR
            kwh=1;
            while kwh==1;
                if kmen2==1 | kmen2==2;
                    prompt={[str1 ' Locally Absent Ring After Year: ']};
                elseif kmen2==3 | kmen2==4;
                    prompt={[str1 ' Locally Absent Ring at Year: ']};
                end;
                def={syrkey};
                dlgTitle=txttit1;
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                yrkey = str2num(answer{1});
                syr =num2str(yrkey);
                
                if kmen2==1 | kmen2==2;
                    if yrkey<yrgo | yrkey>=yrsp;
                        uiwait(msgbox('Must Be an Internal Year','Invalid Entry','modal'));
                        kwh=1;
                    else;
                        kwh=0;
                    end; 
                elseif kmen2==3 | kmen2==4;
                    if yrkey<yrgo | yrkey>yrsp;
                        uiwait(msgbox('Year outside range for series','Invalid Entry','modal'));
                        kwh=1;
                    else;
                        ikey=find(yrx==yrkey);
                        if x(ikey)~=0;
                            uiwait(msgbox([syr ' is not a zero value, pick another year'],'Invalid Entry','modal'));
                        else;
                            kwh=0;
                        end;
                        
                    end; 
                end;
            end; % kwh
            clear kwh
            if kmen2==1; % INSERT AND SHIFT current year and LATER PART FORWARD
                % Store preceding and succeeding parts 
                i1 = yrx<yrkey;
                i2= yrx>=yrkey;
                yrx1=yrx(i1);
                yrx2=yrx(i2)+1;
                x1=x(i1);
                X1=X(i1,:);
                x2=x(i2);
                X2=X(i2,:);
                yrx=[yrx1; yrkey; yrx2];
                x= [x1; 0 ; x2];
                if strcmp(editmode,'Total');
                    X=[X1; 0 ; X2];
                elseif strcmp(editmode,'Partial');
                    X=[X1; 0 0 0; X2];
                end;        
            elseif kmen2==2; % INSERT AND SHIFT EARLIER PART BACK
                % Store preceding and succeeding parts 
                i1 = yrx<=yrkey;
                i2= yrx>yrkey;
                yrx1=yrx(i1)-1;
                yrx2=yrx(i2);
                x1=x(i1);
                X1=X(i1,:);
                x2=x(i2);
                X2=X(i2,:);
                yrx=[yrx1; yrkey; yrx2];
                x= [x1; 0 ; x2];
                if strcmp(editmode,'Total');
                    X=[X1; 0 ; X2];
                elseif strcmp(editmode,'Partial');
                    X=[X1; 0 0 0; X2];
                end;       
            elseif kmen2==3; % REMOVE & KEEP SAME START YEAR
                ikey= find(yrx==yrkey);
                i1=yrx<yrkey;
                i2= yrx>yrkey;
                yrx1=yrx(i1);
                yrx2=yrx(i2)-1;
                x(ikey)=[];
                yrx=[yrx1 ; yrx2];
                X(ikey,:)=[];
            elseif kmen2==4; % REMOVE & KEEP SAME END YEAR
                ikey= find(yrx==yrkey);
                i1=yrx<yrkey;
                i2= yrx>yrkey;
                yrx1=yrx(i1)+1;
                yrx2=yrx(i2);
                x(ikey)=[];
                yrx=[yrx1 ; yrx2];
                X(ikey,:)=[];
            end; % if kmen2==
            clear kmen2 x1 x2 X1 X2 yrx1 yrx2 i1 i2;
            kdone=1;
            
            % Plot series
            subfcn02(x0,yrx0,x,yrx,id{kser},lab1);
        end; % if kdone==1;
        
        
    elseif kmen1==3; % CHANGE REAL RING TO FALSE
        if kdone==1;
            strmsg1='CANNOT STACK CHANGE ON CHANGE: ABORT PREVIOUS CHANGE FIRST';
            uiwait(msgbox(strmsg1,'Message','modal'));
        else;
            yrsp = yrx(end); 
            yrgo = yrx(1);
            syrkey='9999';
            txttit1=['Convert This Real Ring to False'];
            kmen2=menu('Real to False Options',...
                'KEEP SAME START YEAR',...
                'KEEP SAME END YEAR');
            
            % PROMPT FOR YEAR
            kwh=1;
            while kwh==1;
                prompt={'Enter Year'};
                def={syrkey};
                dlgTitle=txttit1;
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                yrkey = str2num(answer{1});
                syr =num2str(yrkey);
                if yrkey<yrgo | yrkey>=yrsp;
                    uiwait(msgbox('Year outside range for series, try again','Invalid Entry','modal'));
                    kwh=1;
                else;
                    kwh=0;
                end;
            end;
            clear kwh
            
            % Add key ring to next 
            ikey = find(yrx==yrkey); % index to year to merge into next 
            x(ikey+1)=x(ikey+1)+x(ikey); % merge ring into following
            X(ikey+1,:)=X(ikey+1,:) + X(ikey,:);
            
            % Remove the key ring
            x(ikey)=[];
            X(ikey,:)=[];
            yrx(ikey)=[];
            nyr = length(yrx);
            
            if kmen2==1; % keep same start year
                yrx = yrgo+((1:nyr)')-1;
                %remeas{kser}=[yrkey]; NO need to remeasure -- just combining existing measurements
            elseif kmen2==2; % keep same end year
                yrx = yrgo+((1:nyr)');
                %remeas{kser}=[yrkey+1];  NO need to remeasure -- just combining existing measurements
            end;
            subfcn02(x0,yrx0,x,yrx,id{kser},lab1);
            kdone=1;
        end;
        
        
    elseif kmen1==4; % CHANGE FALSE RING TO REAL (SPLIT THE RING)
        if kdone==1;
            strmsg1='CANNOT STACK CHANGE ON CHANGE: ABORT PREVIOUS CHANGE FIRST';
            uiwait(msgbox(strmsg1,'Message','modal'));
        else;
            yrsp = yrx(end); 
            yrgo = yrx(1);
            syrkey='9999';
            txttit1=['Convert a false ring in this year to real (split the ring)'];
            kmen2=menu('Real to False Options',...
                'KEEP SAME START YEAR',...
                'KEEP SAME END YEAR');
            
            % PROMPT FOR YEAR
            kwh=1;
            while kwh==1;
                prompt={'Enter Year'};
                def={syrkey};
                dlgTitle=txttit1;
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                yrkey = str2num(answer{1});
                syr =num2str(yrkey);
                if yrkey<yrgo | yrkey>yrsp;
                    uiwait(msgbox('Year outside range for series, try again','Invalid Entry','modal'));
                    kwh=1;
                else;
                    kwh=0;
                end;
            end;
            clear kwh
            
            % Store the key year data and split it in two -- evenly
            ikey=yrx==yrkey;
            xkey=x(ikey);
            xkey=[xkey/2;   (xkey-xkey/2)];
            Xkey=X(ikey,:);
            Xkey = [Xkey/2;  (Xkey-Xkey/2)];
            
            uiwait(msgbox('The ring is automatically split in half; be sure to remeasure to get correct split',...
                'Message','modal'));
            % Rebuild
            if yrkey==yrgo; % key year was first year
                x(1)=[]
                X(1,:)=[];
                x=[xkey; x];
                X=[Xkey; X];
            elseif yrkey==yrsp; % key year was end year
                x(end)=[];
                X(end,:)=[];
                x=[x; xkey];
                X=[X; Xkey];
            else; % key year was internal
                ibefore = find(yrx<yrkey);
                iafter=find(yrx>yrkey);
                x=[x(ibefore); xkey;  x(iafter)];
                X=[X(ibefore,:); Xkey; X(iafter,:)];
            end;
            
            nyr = length(x);
            
            if kmen2==1; % keep same start year
                yrx = yrgo+((1:nyr)')-1;
                remeas{kser}=[yrkey yrkey+1];
            elseif kmen2==2; % keep same end year
                yrx = yrgo+((1:nyr)')-2;
                remeas{kser}=[yrkey-1 yrkey];
            end;
            
            subfcn02(x0,yrx0,x,yrx,id{kser},lab1);
            kdone=1;
        end;
    elseif kmen1==5;% ADD OR REMOVE A STRECH OF MISSING VALUES
        if kdone==1;
            strmsg1='CANNOT STACK CHANGE ON CHANGE: ABORT PREVIOUS CHANGE FIRST';
            uiwait(msgbox(strmsg1,'Message','modal'));
        else;
            yrsp = yrx(end); 
            yrgo = yrx(1);
            syrkey='9999';
            txttit1=['Insert a stretch of missing values beginning after year)'];
            kmen2=menu('Insert missing values option',...
                'KEEP SAME START YEAR',...
                'KEEP SAME END YEAR');
            
            % PROMPT FOR YEAR
            kwh=1;
            while kwh==1;
                prompt={'Enter Year','Enter number of years of missing values'};
                def={syrkey,'1'};
                dlgTitle=txttit1;
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                yrkey = str2num(answer{1});
                nmiss=str2num(answer{2});
                syr =num2str(yrkey);
                if yrkey<yrgo | yrkey>=yrsp;
                    uiwait(msgbox('Year outside range for series, try again','Invalid Entry','modal'));
                    kwh=1;
                else;
                    kwh=0;
                end;
            end;
            clear kwh
            
            % Build NaN block
            xkey =repmat(NaN,nmiss,1);
            Xkey=repmat(NaN,nmiss,size(X,2));
            
            % Rebuild
            ibefore = find(yrx<=yrkey);
            iafter=find(yrx>yrkey);
            x=[x(ibefore); xkey;  x(iafter)];
            X=[X(ibefore,:); Xkey; X(iafter,:)];
            
            nyr = length(x);
            
            if kmen2==1; % keep same start year
                yrx = yrgo+((1:nyr)')-1;
            elseif kmen2==2; % keep same end year
                yrx = flipud(yrsp-((1:nyr)')+1);
            end;
            
            subfcn02(x0,yrx0,x,yrx,id{kser},lab1);
            kdone=1;
        end;
    elseif kmen1==6; % CHANGE A MEASURED VALUE
        
        if kdone==1;
            strmsg1='CANNOT STACK CHANGE ON CHANGE: ABORT PREVIOUS CHANGE FIRST';
            uiwait(msgbox(strmsg1,'Message','modal'));
        else;
            yrsp = yrx(end); 
            yrgo = yrx(1);
            syrkey='9999';
            txttit1=['Change a measured value)'];
            % PROMPT FOR YEAR
            kwh=1;
            while kwh==1;
                prompt={'Enter Year'};
                def={syrkey};
                dlgTitle=txttit1;
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                yrkey = str2num(answer{1});
                syr =num2str(yrkey);
                if yrkey<yrgo | yrkey>yrsp;
                    uiwait(msgbox('Year outside range for series, try again','Invalid Entry','modal'));
                    kwh=1;
                else;
                    kwh=0;
                end;
            end;
            clear kwh
            ikey=yrx==yrkey;
            xkey=x(ikey);
            Xkey=X(ikey,:);
            % Prompt for value
            prompt={'Enter Value (mm)'};
            def={num2str(xkey)};
            dlgTitle=txttit1;
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            xreplace= str2num(answer{1});
            
            % Make same split as befor on partial
            if size(X,2)==3;
                p1 = Xkey(1)/Xkey(3);
                p2= 1-p1;
                Xreplace=[p1*xreplace p2*xreplace xreplace];
            else;
                Xreplace=xreplace;
            end;
            
            % Replace
            x(ikey)=xreplace;
            X(ikey,:)=Xreplace;
            remeas{kser}=yrkey;
            
            
            subfcn02(x0,yrx0,x,yrx,id{kser},lab1);
            kdone=1;
        end;
        
    elseif kmen1==7; % ABORT THE CHANGE AND START OVER
        yrx=yrx0;
        x=x0;
        X=X0;
        remeas{kser}=[];
        dspan{kser}=V.span{kser};
        kdone=0;
        close(1);
        figure(1);
        subfcn01(x0,yrx0,id{kser},lab1);
    elseif kmen1==8; % ACCEPT THE CHANGE AND QUIT
        dspan{kser}=[yrx(1) yrx(end)];
        V.span = dspan;
        V.remeasure=remeas;
        d{kser}=[yrx X];
        V.data=d;
        
        kwh1=0;
    end; % if kmen1==
end; % while


function subfcn01(x0,yrx0,id,lab1);
% Plot series
strrng=[num2str(yrx0(1)) '-' num2str(yrx0(end))];
hp2=plot(yrx0,x0);
title([id ':  before editing (' strrng ')']);
xlabel('Year');
ylabel(lab1); 
grid;
zoom xon;


function subfcn02(x0,yrx0,x,yrx,id,lab1);
% Plot series
strrng0=[num2str(yrx0(1)) '-' num2str(yrx0(end))];
strrng1=[num2str(yrx(1)) '-' num2str(yrx(end))];
hp2=plot(yrx0,x0,yrx,x);
set(hp2(2),'LineWidth',2);
title([id ' before and after editing']);
xlabel('Year');
ylabel(lab1); 
legend(['Before (' strrng0 ')'],['After (' strrng1 ')']);
grid;
zoom xon;

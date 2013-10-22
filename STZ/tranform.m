function tranform
% tranform;
% tranform: interactive power transformation of ring-width series
% Last revised 1-3-01
%
% One of sequence of functions for tree-ring standardization.  Typically, follows appliction of rwlinp.m and precedes 
% crvfit.m.  Applies function sprdloc.m to log-transform or power transform ring widths. Details of transformation
% method and options in comments of sprdloc.m.
%
%*** IN ****************************************
%
% User prompted for name of .mat file storing the ringwidth data and associated years and names.  Ringwidth data 
% previously put in this file with rwlinp.m
%
% Critical input data is:
% X (? x 1)r vector with measured ring widths
% nms (? x 
% yrs
%
%*** OUT ***********************************
%
% User prompted for name of .mat file storing transformation information. Can store in a new file, or append to the 
% input .mat file.  Typically append.  
%
% Output variables produced or modified:
%  A {? x 8} summary of transformations for ? cores
%   {1} c (1 x 1)r -- shift increment (quantity x+c) is transformed
%   {2} p(1 x 1)r -- power of transformation 
%       p==1 no transform
%       p==0 log transform
%       p otherwise: (x+c) raised to power p
%   {3} a (1 x 1)r first coefficient for matching  ( y = a + bT(x+c))
%   {4} b (1 x 1)r second coefficient for matching 
%   {5} lenseg (1 x 1)r size of batches for spread-vs-level computations
%   {6} shift (1 x 1)r of batches for spread-vs-level computations
%   {7} khow (1 x 1)r -- how transformation selected)
%       ==1 automatic choice of log-transform (base 10)
%       ==2 automatic choice of commonly used power compatible with slope of spread vs level (see notes for sprdloc.m)
%       ==3 automatic choice of power computed from slope of spread vs level
%       ==4 interactive fitting, but accepting power computed from slope of spread-vs-level
%       ==5 interacitve fitting, specifying some other power
%   {8} eqnstr (1 x ?)s equation for transformation (without a,b,c parameters built in)
%   tstatus (? x 1)i transformation status of individual series, where ? is number of series
%       ==0 never transformed
%       ==1 transformed, but not revised (see notes)
%       ==2 revised transformation (see note)
% Fwhen {? x 4} history of chronology development (date, hour, minute that various programs rin
%   {1,:} rwlinp
%   {2,:} tranform
%   {3,:} crvfit
%   {4,:} corei
%   {5,:} treei
%   {6,:} sitei
%
%   The four cells for each row hold
%       function : e.g., "rwlinp"
%       date&time: e.g., Jan 01 01  10 3 == Jan 1 2001, hour 10, minute 3
%       input file:  e.g., "ste.rwl"
%       output file:  e.g., "ste.mat"
%  
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED (FROM C:\MLB\STZ UNLESS NOTED)
%*** TOOLBOXES NEEDED
%
%*** NOTES
%
% tstatus.  Facilitates revising transformations.  Say all series have been transformed in one of the automatic modes. User 
% then might want to sequentially view and modify the individual transformations.  The automatic run would set tstatus=1 for 
% all series.  As each is modified, the status changes to 2.  When all have been modified, tstatus for all is changed 
% from 1 to 2.  
%
% Local re-transformation.  This is series by series manual interactive revising of transformations.  Transformations can 
% be individualized.  But a single setting for window length and shift (legseg and mlap) applies to all series.
%
% Window size and shift for call to sprdloc.  The parameters lenseg and mlap are prompted for.  Defaults are 30 yr and 
% 10 yr.  For very short series, these defaults might be too large.  If the series length is less than or equal to 3 times 
% lenseg, alternative, smaller, values for lenseg and mlap are computed for those series. Roughly, lenseg is set to 1/3 the 
% series length, and mlap to 1/3 of lenseg.


% Close any open windows
close all

% Prompt for name of .mat file with ring widths, core ids, and 
% year information
[flmat,path1]=uigetfile('*.mat','Input .mat ringwidth storage file');
pf1=[path1 flmat];
flold=flmat;

% Load the .mat storage file
eval(['load ' pf1]);

% The .mat storage file should contain X, nms, and yrs
% Also may contain growth trend data  G, Gnms, Gyrs
if ~(exist('X')==1) | ~(exist('nms')==1) | ~(exist('yrs')==1),
  error('The selected .mat file does not contain X, nms and yrs');
end

% Get the number of cores, which equals the row size of  matrix
% of core ids
ns=size(nms,1);

% Allocate history
Lwhen=exist('Fwhen','var');
if ~Lwhen;
    Fwhen = cell(8,4);
end;

% Store history information 
Fwhen{2,1}='tranform';
Fwhen{2,3}=flmat;
clear flmat



%***********************************************************************

% Prompt for batch parameters for call to sprdloc
prompt={'Enter window length:','Enter the offset:'};
def={'30','5'};
dlgTitle='Sliding window for defining batches';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
lenseg=str2num(answer{1});
mlap=str2num(answer{2});


% Prompt for allowing use of subperiod for choosing transformation
ksubpd = questdlg('Allow point/click on subperiod for picking transformations?');
switch ksubpd;
case 'Yes';
    Period='Manual';
otherwise;
    Period='Auto';  % Fit period specified in S(10:12).  If no fit has been specified yet, uses full period
end;


% Find out if any previous transformations

if exist('A','var')~=1; % If no previous transformation info initialize storage 
   A = cell(ns,8);
   tstatus = zeros(ns,1);
   Lnew=1;
else;
    switch Period;
    case 'Auto';
        
        Lnew=0;
        quest= questdlg('Globally re-transform series?');
        switch quest;
        case 'Cancel';
            close all;
            error('Pick Yes or No for quest on Globally re-transform');
        case 'Yes';
            Lglobal=1;
        case 'No';
            Lglobal=0;
        end;
    case 'Manual';
        Lglobal=0;
    otherwise;
    end; % switch Period;
    
end;


txtin{2}='Ring width';
txtin{3}= 'mm x 100';


if Lnew==1 | Lglobal==1;  % Handle first transformation or global re-transformation
    fprev=[];
    % Prompt for automatic transformation option
    kmen3=menu('Choose method for transforming all series',...
        'Log10',...
        'Conservative common power',...
        'Power compute from slope of spread vs level');
    if kmen3==1; % Log10
        kopt=1; % call option for srdloc
    elseif kmen3==2;
        kopt=2;
    elseif kmen3==3;
        kopt=3;
    end;
    
    for n =1:ns; % loop over cores
        
        % Get the ring widths
        if exist('S','var')==1 & S(n,12)~=0;
            i1=S(n,12);
            yrgo=S(n,10);
            yrsp=S(n,11);
        else;
            i1=yrs(n,3);
            yrgo = yrs(n,1);
            yrsp = yrs(n,2);
        end;
        nyr = yrsp-yrgo+1; 
        
        if nyr<3*lenseg;
            lenseg0 = floor(nyr/3);
            mlap0=floor(lenseg0/3);
        else;
            lenseg0=lenseg;
            mlap0=mlap;
        end
        
        i2= i1 + nyr -1;
        x  = X(i1:i2);
        yrx=(yrgo:yrsp)';
        
        % Store the text input required by sprdloc.m
        txtin{1}=nms(n,:); % id
        
               
        % Transform
        D=[yrx x];
        [Y,f]=sprdloc(D,lenseg0,mlap0,txtin,fprev,kopt);
        A{n,1}=f{1};
        A{n,2}=f{2};
        A{n,3}=f{3};
        A{n,4}=f{4};
        A{n,5}=lenseg0;
        A{n,6}=mlap0;
        A{n,7}=f{6};
        A{n,8}=f{5};
    end;
    tstatus(:)=1; % transformation status
else;  % Not first time this data has been transformed
    
    % OPTIONAL LOCAL TRANSFORMATION
    
    % If all series previously re-transformed re-set flag to 1 for all
    if all(tstatus==2);
        tstatus(:)=1;
    end
    
    kopt(1)=4; % manual option for call to sprdloc
    
    % Build string matrix of status into names
    ts =cellstr([nms repmat('-',ns,1)  num2str(tstatus)]);
    ts{ns+1}='Quit';
    knew=1;
    
    kwh1=1; % while working on the individual series
    while kwh1;
        
        if knew==1;
            n=menu('Choose ',ts);
        end;
        
        % Break out of kwh1 loop if chose quit
        if n==(ns+1);
            break;
        end;
                
        % Get and store previous transformation parameters
        fprev=cell(1,6);
        fprev{1}=A{n,1}; % c
        fprev{2}=A{n,2}; % p
        fprev{3}=A{n,3}; % a
        fprev{4}=A{n,4}; % b
        fprev{5}=char(A{n,8}); % eqnstr    
        fprev{6}=A{n,7}; % khow
        
        % Get the ring widths
        if exist('S','var')==1 & S(n,12)~=0;
            i1=S(n,12);
            yrgo=S(n,10);
            yrsp=S(n,11);
        else;
            i1=yrs(n,3);
            yrgo = yrs(n,1);
            yrsp = yrs(n,2);
        end;
        nyr = yrsp-yrgo+1; 
        
        if nyr<=3*lenseg;
            lenseg0 = floor(nyr/3);
            mlap0=floor(lenseg0/3);
        else;
            lenseg0=lenseg;
            mlap0=mlap;
        end
        
        i2= i1 + nyr -1;
        x  = X(i1:i2);
        yrx=(yrgo:yrsp)';
        
        % Full-length variables
        i1orig=yrs(n,3);
        yrgoorig = yrs(n,1);
        yrsporig = yrs(n,2);
        nyrorig = yrsporig-yrgoorig+1; 
        yrxorig=(yrgoorig:yrsporig)';
        i2orig=i1orig+nyrorig-1;
        xorig=X(i1orig:i2orig);
        
        % Store the text input required by sprdloc.m
        txtin{1}=nms(n,:); % id
        
        %------- Optional use of subperiod for selecting transformation.
        switch Period;
        case 'Manual';
            figure(3); % time series
            hp2=plot(yrx,x,yrxorig,xorig,':');
            ylabel([txtin{2} '(' txtin{3} ')']);
            xlabel('Year');
            legend('Stored subperiod','Full-period');
            
            title(txtin{1});
            
            kmen1=menu('Choose period for fitting transformation',...
                'Full period',...
                'Click on ends',...
                'Stored subperiod (solid blue line)');
            if kmen1==1; % full period
                xx=xorig;
                yrxx=yrxorig;
            elseif kmen1==2; % click on ends
                
                jh1=msgbox('Click at two corner points of fit period',' ');
                pause(1);
                close(jh1);
                figure(3);
                [segx,segy]=ginput(2);
                [erflg,xind1,xind2]=erchk(yrxorig,min(segx),max(segx));
                
                yrxx=yrxorig(xind1:xind2); % year vector for selected fit interval
                xx=xorig(xind1:xind2);   % ring-width data for the selected fit interval
            elseif kmen1==3; % use subperiod stores in S(i,10:12);
                xx=x;
                yrxx=yrx;
            end;
        case 'Auto'; % use default period to select transform. This is full period if S does not exist or period
            xx=x;
            yrxx=yrx;
          
        otherwise;
        end; % switch Period;
        
        
        
        % Transform
        D=[yrxx xx];
        [Y,f]=sprdloc(D,lenseg0,mlap0,txtin,fprev,kopt);
        
        close all;
        
        figure(1); % time series
        xmed=median(D(:,2));
        hp1=plot(yrxx,D(:,2),yrxx,Y(:,2),[min(yrxx) max(yrxx)],[xmed xmed]);
        ylabel([txtin{2} '(' txtin{2} ')']);
        xlabel('Year');
        legend('Oringinal RW','Transformed RW','Median');
        strtit=strtok(ts{n},'-');
        title([strtit ';  ' f{5}]);
        
        figure(2) ; % box
        G=[D(:,2) Y(:,2)];
        boxplot(G);
        set(gca,'XTickLabel',{'Original','Transformed'});
        title('Effect of Transformation on Distribution');
        ylabel([txtin{2} '(' txtin{2} ')']);
        xlabel('Variables');
        
        figure(1);
        
       
        
        kmen4=menu('Choose','Accept the transformation','Try again');
        if kmen4==1;
            % Store the transformation info
            A{n,1}=f{1};
            A{n,2}=f{2};
            A{n,3}=f{3};
            A{n,4}=f{4};
            A{n,5}=lenseg0;
            A{n,6}=mlap0;
            A{n,7}=f{6};
            A{n,8}=f{5};
            % Update status and new series flag
            tstatus(n)=2;
            knew = 1;
            % Revise the menu text
            s1=ts{n}
            lens1=length(s1);
            s1(lens1)='2';
            ts{n}=s1;
        else;
            knew=0;
        end;
        
    end; % while kwh1
end; % 


c=clock;
c=num2str(c(4:5));
d=date;
Fwhen{2,2}=[d ', ' c];
% Save the vectors in a .mat file
nsv=menu('Save variables?','Add to original .mat storage file?',...
   'Store X, yrs, nms, A, ,tstatus, Fwhen in a new file?','QUIT');
if nsv==1; % add variables to the original input .mat file
    Fwhen{2,4}=flold;   
   eval(['save ' pf1 ' A Fwhen tstatus -append']);
elseif nsv==2,
   [ofmat,path2]=uiputfile('*.mat','new .MAT file to store X,nms,yrs,A,Fwhen, tstatus: ');
   pf2=[path2 ofmat];
   Fwhen{2,4}=ofmat;   
   eval(['save ' pf2 ' X yrs nms A Fwhen tstatus ']);
end
% End of file

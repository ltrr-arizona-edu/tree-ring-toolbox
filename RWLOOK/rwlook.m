function rwlook
% rwlook: bi-variate graphical and statistical quality-control of crossdating
% CALL: rwlook;
%
% Last revised 10-21-02
%
%********** INPUT
%
% No input arguments
% User prompted to pick two series for comparison. 
%
%********* NOTES
%
% First series is the test series.  Second series is the reference series.
% Several modes of use are possible. You could compare:
%  a questionably dated .rw file with a well-dated .rw file
%  a questionably dated .rw file with a master chronology
%  a tentatively dated chronology with a nearby chronology
%
% rwlook has some editing capability.  Only the test series can be edited.  You
% can check by trial and error the effect of adding locally absent rings, changing
% true ring to false ring, etc.
%
% If .rw files are to be used in running rwloook, you need to
% run rwfiles.m first to build and store a cell variable with core ids. 
% rwlook uses the stored cell variable to make its initial screen menu of .rw files.
% If user revises .rw data in any way to change the set of .rw filenames that should
% appear in the screen menu, user must run rwfiles.m to update the cell variable before
% running rwlook.m
%
% Note that prefixes.mat is hardcoded as the .mat file that holds the cell variable 
% produced by rwfiles.m (see "HARD CODED" below)
%
%
% Display of the raw time series (ring widths or indices).  Time series plots of the 
% test series and reference series on the same axes are part of the graphical output. 
% If the series are of different units or scale (e.g., indices and ring widths), some
% some adjustments are needed to ensure that the plots can be readily compared.  Care 
% also has to be taken that smoothed values of the adjusted series will not approach 
% zero or be less than zeros, because in computing the "ring-width change" series, the 
% smoothed series is used as a denominator in high-pass filtering. Thus the following 
% adjustment:
%   If both series are from .rw files, the plots are the ring widths as measured, 
%      in hundredths of mm
%   If the test series is from a .rw file and the the other reference series from a 
%      .dat file, the reference series is scaled and shifted as needed to make the 
%      plots more readable.  The reference series is scaled to the same variance as 
%      the test series, and is shifted to have the same mean as the test series.  In 
%      doing this, it is possible that the resulting adjusted reference series will 
%      have some values below zero.  These values would wreak havoc in later filtering 
%      by possible leading to division by zero or by some some value near zero.  To 
%      handle this problem, if any of the adjusted values of the reference series are
%      below zero, the adjusted series is shifted vertically so that its minimum value 
%      is the same as that of the test series.
%   If both series are from .dat files, the series are plotted as input, except that
%      any negagive values are set to zero.  This fix is needed because negative indices 
%      sometimes crop up in ARMA-residual chronologies
% ******************   U-W FUNCTIONS NEEDED  ***************
%
%----- from \rwlook\
% rwfiles.m -- needed in pre-processing
% rwread.m
% rwchng.m
% rwedit.m
% rwwrite.m
%  
%------- from \mlb\
% pullseg1.m
% kendtau.m
% signtest.m
% hypernew.m
% quantile.m
% pltext.m
%
% signtbl.mat
% cona12.mat   -- Conover table A12 for testing tau statistic
%
%
%***********************************************************


clear;
close all; 

%-- Message box to define series 1 and 2 to user

hmsg1={'In the terminology of this function,'...,
        'Series 1 is the test series, and ',...
        'Series 2 is the reference series',...
        ' ',...
        'The test series is the series whose dating or measurements might be in',...
        'doubt.  The dating of the reference series should be unequivocal, and its',...
        'measurements accurate.'};
uiwait(msgbox(hmsg1,'Message','modal'));

kwhsource=1;
while kwhsource==1;
    ksource=menu('Choose one',...
        'Series #1 and #2 from same rwmeas storage file',...
        'Series #1 and #2 from different rwmeas storage files',...
        'Series #1 from rwmeas storage file, Series #2 from ITRDB -format .crn file',...
        'Series #1 and Series #2 from separate .crn files',...
        'Series #1 and Series #2 from .rw files',...
        'Series #1 from .rw file, Series #2 from .crn file',...
        'Series #1 and #2 from separate 2-column ascii files');
    if ~any(ksource==[1 2 3 4]);
        uiwait(msgbox({'Only first three menu choices supported so far','Message'},'modal'));
    else;
        kwhsource=0;
    end;
end;


%---- Option for skipping loading of prefixes.mat

if any (ksource==[5 6]);
    kpre=questdlg('Do you need to load prefixes.mat');
    if strcmp(kpre,'Yes');
        uiwait(msgbox({'You are assumed to have run rwfiles.m on a group of . rw files.',...
                ' That would have produced the file prefixes .mat in the current working directory','Message'},'modal'));
        load prefixes; % cell array of core prefixes in pre
    else;
    end;
else;
    kpre='No';  % will not need to load prefixes.mat
end;
if any (ksource==[1 2 3]);
    kmeas='Yes'; % one or both from an rwmeas storage file
else;
    kmeas='No';
end;


%-- OPTIONAL LOADING OF .MAT STORAGE FILE

if strcmp(kmeas,'Yes');
    if any(ksource==[1 3]); % either both series from same an rwmeas storage, or only first series from and rwmeas storage
        [file3,path3]=uigetfile('*.mat','Input .mat storage file with rwmeas-formatted ring-width data');
        pf3=[path3 file3];
        eval(['load ' pf3 ' XT XEL;']);
    elseif  ksource==2; % need two rwmeas storage files
        [file3,path3]=uigetfile('*.mat','Input .mat storage file with first series');
        pf3=[path3 file3];
        eval(['load ' pf3 ' XT XEL;']);
        XT1=XT;
        XEL1=XEL;
        [file3a,path3a]=uigetfile('*.mat','Input  .mat storage file with second series');
        pf3a=[path3a file3a];
        eval(['load ' pf3a ' XT XEL;']);
        XT2=XT;
        XEL2=XEL;
    else;
    end;
end;



% Initialize matrix to store log information on editing
E=blanks(15);

% Load lookup table for signif test of Kendall's tau
if ~exist('cona12');
	load cona12;
end


%-- PROMPT FOR SOURCES OF SERIES 1 AND 2

if strcmp(kmeas,'Yes');
    kwh5=1;
    while kwh5;
        if ksource==1; % series 1 and 2 from same rwmeas storage file
            type1='meas'; type2='meas';
            kmode=menu('Choose mode of run',...
                'Compare EWW of two cores',...
                'Compare LWW of two cores',...
                'Compare measured TWW of two cores',...
                'Compare computed TWW of two cores',...
                'Compare computed TWW with measured TWW of same core');
            if kmode==1 | kmode==2 | kmode==4;
                XX=XEL;
                YY=XEL;
                kwh5=0;
            elseif kmode==3;
                XX=XT;
                YY=XT;
                kwh5=0;
            elseif kmode==5;
                uiwait(msgbox('Make another choice: rwcomp.m is the function for this!','Message','modal'));
            end;
            
            if length(XX.id)<2 ;
                uiwait(msgbox('Only one series of desired type in the rwmeas storage file','Message','modal'));
                kkill=questdlg('Abort the run');
                if strcmp(kkill,'Yes');
                    error('Aborted');
                else;
                end;
            end;
            
            
             if kmode==1; 
                icol =2;
            elseif kmode==2;
                icol=3;
            elseif kmode==3;
                icol=2;
            elseif kmode==4;
                icol=4;
            end;
        elseif ksource==3;  % series 1 from rwmeas, series 2 from crn
            kwh5=0;
            type1='meas'; type2='crn';
            kmode=menu('Choose mode of run',...
                'EWW',...
                'LWW',...
                'Measured TWW',...
                'Computed TWW');
            if kmode==1 | kmode==2 | kmode==4;
                XX=XEL;
            else;
                XX=XT;
            end;
            % Get .crn data
            [file4,path4]=uigetfile('*.crn','Input .crn file with reference chronology');
            pf4=[path4 file4];
            [Y1orig,snull,yr]=crn2vec2(pf4);
            Y1orig=[yr Y1orig];
            Yfnorig=file4;
            YY.id{1}=Yfnorig;
             if kmode==1; 
                icol =2;
            elseif kmode==2;
                icol=3;
            elseif kmode==3;
                icol=2;
            elseif kmode==4;
                icol=4;
            end;
            clear snull yr;
        elseif ksource==2; % series 1 and 2 from separate rwmeas storage files
            kwh5=0;
            type1='meas'; type2='meas';
            kmode=menu('Choose mode of run',...
                'EWW',...
                'LWW',...
                'Measured TWW',...
                'Computed TWW');
            if kmode==1 | kmode==2 | kmode==4;
                XX=XEL1;
                YY=XEL2;
            else;
                XX=XT1;
                YY=XT2
            end;
            if kmode==1; 
                icol =2;
            elseif kmode==2;
                icol=3;
            elseif kmode==3;
                icol=2;
            elseif kmode==4;
                icol=4;
            end;
        
        end; % ksource==
    end; % kwh5
else;  % Neither series if from an rmeas storage file
    if ksource==4; % both from .crn files
        type1='crn'; type2='crn';
        % Get .crn data, first series
        [file4a,path4a]=uigetfile('*.crn','Input .crn file with test chronology');
        pf4a=[path4a file4a];
        [X1orig,snull,yr]=crn2vec2(pf4a);
        X1orig=[yr X1orig];
        Xfnorig=file4a;
        XX.id{1}=Xfnorig;
        % Get .crn data, second series
        [file4,path4]=uigetfile('*.crn','Input .crn file with reference chronology');
        pf4=[path4 file4];
        [Y1orig,snull,yr]=crn2vec2(pf4);
        Y1orig=[yr Y1orig];
        Yfnorig=file4;
        YY.id{1}=Yfnorig;
        clear snull yr;
    elseif ksource==5;
        type1='rw'; type2='rw';
    elseif ksource==6;
        type1='rw'; type2='dat';
    elseif ksource==7;
        type1='dat'; type2='dat';
    end;
end;
clear kmen1;

% Check that will have at least one series in list to pick series 1 from and one series in list to pick series 2 from
if any(ksource==[1 2 3])    & (isempty(XX.id{1}) | isempty(YY.id{1}));
    error('No series of the desired type (e.g., EWW) in one or both the rwmeas storage files');
end;


kk2=1;
while kk2==1;
    
    if  ~strcmp(kmeas,'Yes'); % neither series from rwmeas formatted .mat file
        
        %--- Read in first series -- from a .rw file
        if strcmp(type1,'rw');
            [X1,perx,whenx,Xfn]= rwread(rwpre,'first');
            % X1 is now a 2-col vector with year in column 1
        elseif strcmp(type1,'dat'); % first series is from a .dat file
            [file1,path1]=uigetfile('*.dat','Input file with series #1');
            pf1=[path1 file1];
            eval(['load ' pf1 ' -ascii;']);
            eval(['X1= ' strtok(file1,'.') ';']); 
            Xfn = file1;
            perx=[]; whenx=[];
         elseif strcmp(type1,'crn'); % first series is from a .crn file
            X1=X1orig; 
            Xfn = Xfnorig;
            perx=[]; whenx=[];
        end;
        
        %-- Read in second series
        if strcmp(type2,'rw'); % reference series from a .rw file
            if strcmp(type1,'dat');
                error('Not allowed to have series 1 as .dat and series 2 as .rw');
            end
            [Y1,pery,wheny,Yfn]=rwread(rwpre,'second');
        elseif strcmp(type2,'dat'); % from a 2-column .dat file
            [file2,path2]=uigetfile('*.dat','Input file with series #2');
            pf2=[path2 file2];
            eval(['load ' pf2 ' -ascii;']);
            eval(['Y1= ' strtok(file2,'.') ';']); 
            Yfn = file2;
            pery=[]; wheny=[];
        elseif strcmp(type2,'crn'); % second series is from a .crn file
            Y1=Y1orig; 
            Yfn = Yfnorig;
            pery=[]; wheny=[];
        end;
        
    else; % test series (and maybe ref series too, from rwmeas file)
        if kmode~=5;  % Not comparing computed with measured total width of same core
            C1=XX.id';
            nlen = length(C1);
            Lin = logical(zeros(nlen,1));
            C2=cellstr(repmat('-N',nlen,1));
            nallow=[1 1];
            strmenu='Toggle to Y to pick first series';
            [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu);
            C1{Lout}=['**' C1{Lout} '**'];
        else; % comparing computed with measured total width for same core
            % Must identify cores for which have both partial and total width measurements
                   
        end;
        
        if any(kmode==[1 2 4]) % EWW, LWW or computed total width of two cores
            X1=[XX.data{find(Lout)}(:,1)    XX.data{find(Lout)}(:,icol)];
            perx=[];
            whenx=[];
            Xfn=XX.id{find(Lout)};
            
            if strcmp(type2,'meas');
                strmenu='Toggle to Y to pick second series';
                nallow=[1 1];
                pery=[];
                wheny=[];
                if ksource==2; % series 2 from differenc rwmeas file
                    C1b=YY.id';
                    nlenb = length(C1b);
                    Linb = logical(zeros(nlenb,1));
                    C2b=cellstr(repmat('-N',nlenb,1));
                    [C3b,Loutb] = menucell (C1b,C2b,Linb,nallow,strmenu);
                    C1b{Loutb}=['**' C1b{Loutb} '**'];
                    Y1=[YY.data{find(Loutb)}(:,1)    YY.data{find(Loutb)}(:,icol)];
                    Yfn=YY.id{find(Loutb)};
                else;
                    Lin = logical(zeros(nlen,1));
                    C2=cellstr(repmat('-N',nlen,1));
                    [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu);
                    Y1=[YY.data{find(Lout)}(:,1)    XX.data{find(Lout)}(:,icol)];
                    Yfn=YY.id{find(Lout)};
                end;
                
            elseif strcmp(type2,'crn');
                Y1=Y1orig;
                Yfn=Yfnorig;
                pery=[];
                wheny=[];
            else;
                error(['type2 == ' type2 ' invalid']);
            end;
        
        elseif kmode==3; % measured TWW of two cores
             X1=[XX.data{find(Lout)}(:,1)    XX.data{find(Lout)}(:,2)];
            perx=[];
            whenx=[];
            Xfn=XX.id{find(Lout)};
            
            if strcmp(type2,'meas');
                strmenu='Toggle to Y to pick second series';
                Lin = logical(zeros(nlen,1));
                C2=cellstr(repmat('-N',nlen,1));
                nallow=[1 1];
                [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu);
                Y1=[YY.data{find(Lout)}(:,1)    XX.data{find(Lout)}(:,2)];
                pery=[];
                wheny=[];
                Yfn=YY.id{find(Lout)};
            elseif strcmp(type2,'crn');
                Y1=Y1orig;
                Yfn=Yfnorig;
                pery=[];
                wheny=[];
            else;
                error(['type2 == ' type2 ' invalid']);
            end;
            
                   
        elseif kmode==5; % computed TWW with measured TWW of same cores
            X1=[XX.data{find(Lout)}(:,1)    XX.data{find(Lout)}(:,4)]; % computed total width
            perx=[];
            whenx=[];
            Xfn=XX.id{find(Lout)};
            error('NOT YET CODED: Need code to compare selected ID of XEL with ids of XT and get same');
        end;
        
        
    end;
    
    %-- If reference series from a .dat file and test series from a .rw file, 
    %-- adjust reference series to same mean and standard dev as test series
    if (strcmp(type2,'dat') | strcmp(type2,'crn'))     & (strcmp(type1,'rw') |  strcmp(type1,'meas')  );
        s1 = std(X1(:,2));
        mn1 = mean(X1(:,2));
        ytemp =zscore(Y1(:,2));
        Y1(:,2)=mn1 + ytemp * s1;
        clear s1 mn1 ytemp;
        
        %--- If any negative values in the scaled reference series,
        % adjust level of reference series so that its lowest value is equal to the
        % lowest value of the floating series
        Lneg = any(Y1(:,2)<0);
        if Lneg;
            xlow = min(X1(:,2));
            ylow = min(Y1(:,2));
            yraise = xlow-ylow;
            Y1(:,2)=Y1(:,2)+ yraise;
        end;
        clear yraise xlow ylow Lneg;
    end;
    
    
    % If test series and reference series are both from .dat files, I assume 
    % both are index chronologies.  In this case neither series is scaled, so that
    % the "ring-width" plots are the chronology indices. But to avoid possible problems 
    % with division by near-to-zero, any negative indices are changed to zero. Such
    % negative indices sometimes crop up in "residual" chronologies.
    if (strcmp(type1,'dat') & strcmp(type2,'dat')) (strcmp(type1,'crn') & strcmp(type2,'crn'));
        Lneg = X1(:,2)<0;
        if any(Lneg);
            X1(Lneg,2)=0;
        end;
        Lneg = Y1(:,2)<0;
        if any(Lneg);
            Y1(Lneg,2)=0;
        end;
    end;
    
    % For convenience, set X to the test series, Y to the reference series
    X=X1;
    Y=Y1;
    
    
    %--- Truncate file names for ease of reading on plots
    
    % Questionable series-X
    if strcmp(type1,'rw');
        n1=length(Xfn);
        Xfnlab=Xfn;
        Xfnlab(n1-2:n1)=[];  % get rid of '.rw'
        cn=':';
        cn=cn(:,ones(length(Xfnlab),1));
        cni=find(cn==Xfnlab);
        if cni>=0,
            Xfnlab(1:cni)=[]; % gets rid of leading drive info, like "B:"
        end;
    elseif strcmp(type1,'dat') | strcmp(type1,'crn');
        Xfnlab = strtok(Xfn,'.');
    elseif strcmp(type1,'meas');
        Xfnlab=Xfn;
    end;
    
    
    % Reference Series-- Y
    if strcmp(type2,'rw');
        Yfnlab=Yfn;
        n1=length(Yfnlab);
        Yfnlab(n1-2:n1)=[];  % get rid of '.rw'
        cn=':';
        cn=cn(:,ones(length(Yfnlab),1));
        cni=find(cn==Yfnlab);
        if cni>=0,
            Yfnlab(1:cni)=[]; % gets rid of leading drive info, like "B:"
        end;
    elseif strcmp(type2,'dat')  | strcmp(type2,'crn');
        Yfnlab=strtok(Yfn,'.');
     elseif strcmp(type2,'meas');
        Yfnlab=Yfn;    
    end;
    
    
    kfirst=1;
    
    kk1=1; % while control for level-0 menu
    while kk1<5; %  kk1=5 would exit program
        
        if kfirst==0; % Not first time thru, so want option of things to do
            kk1=menu('Choose one',...
                'Analyze ',...
                'Edit rw ',...
                'Overwrite old .rw with revised ringwidths',...
                'Make edit window',...
                'Finished working with this pair of series',...
                'Quit program');
        else; % kfirst not zero; first pass thru; will want to analyze
            kk1=1;
        end
        
        if kk1==1;  % analyze
            rwlist(X,Y);
            % Set vectors with rw and years
            x=X(:,2);  % series
            yrx=X(:,1);  % year vector
            y=Y(:,2);  % series 2
            yry=Y(:,1);  % year vector
            nx=length(x);
            ny=length(y);
            
            % Transform ring width
            w=rwchng(x);
            z=rwchng(y);
            yrw=yrx(2:nx);
            yrz=yry(2:ny);
            
            k1=1;        % while control for level-1 menu
            while k1~=4;  % if k1==4, will return to main menu
                if kfirst==0;
                    k1=menu('Choose one: ','tau',...
                        'Zoom RW','Zoom D','Return to main menu');
                elseif kfirst==1;
                    k1=1;
                end
                
                
                if k1==1;   % Initial sign and kendal-tau tests and plots
                    % do this only on transformed data
                    
                    if kfirst==0;
                        prompt={'Enter the Segment Length (yr):','Enter the offset(yr):'};
                        def={'50','5'};
                        %def={'30','1'};
                        tit1='Tau analysis settings';
                        lineNo=1;
                        answer=inputdlg(prompt,tit1,lineNo,def);
                        len=str2num(answer{1});
                        offset=str2num(answer{2});
                        
                    else
                        % here
                        len=50;
                        offset=5;
                    end
                    
                    txttau1 = ['\itM\rm=' int2str(len) ' yr,  '];
                    txttau2 = ['\Deltat = ' int2str(offset) ' yr'];
                    txttau=[txttau1 txttau2];
                    
                    % Find overlap between transformed rw series
                    yrso=[max([yrw(1) yrz(1)]) min([yrw(nx-1) yrz(ny-1)])];
                    
                    
                    
                    yro=(yrso(1):yrso(length(yrso)))';
                    L3w=yrw>=yrso(1) & yrw <=yrso(length(yrso));
                    L3z=yrz>=yrso(1) & yrz <= yrso(length(yrso));
                    
                    [yreh,H]=pullseg1(w(L3w),yrw(L3w),len,offset); % for w
                    [yreb,B]=pullseg1(z(L3z),yrz(L3z),len,offset); % for z
                    
                    [mH,nH]=size(H);
                    v=zeros(nH,1);
                    
                    for i=1:nH;
                        h=H(:,i);
                        b=B(:,i);
                        [T,tau]=kendtau(h,b);
                        v(i)=tau;
                        
                    end
                    vv=v;
                    
                    % 0.95 and 0.99 one-sided confidence levels for tau
                    % Revision 1-31-94
                    % One-sided test appropriate.  H0 is that tau=0.  H1 is that
                    %   tau>0.  Function not written to test for significance of
                    %   negative correlation, because unlikely to be of interest in
                    %   comparing ring-width series.
                    % Tcrit=table1(cona12,len); % get a row from the lookup table -- OBSOLETE
                    Tcrit=interp1(cona12(:,1),cona12(:,2:size(cona12,2)),len); % REVISED 10-23-02
                    tau95=Tcrit(2)/(len*(len-1)/2);  % for 95% signif level;
                    % Denominator converts T to tau
                    tau99=Tcrit(4)/(len*(len-1)/2);  % 99% signif level
                    
                    
                    figure(1);
                    plot(yreh,v,yreh,tau95(ones(nH,1),:),'m--',...
                        yreh,tau99(ones(nH,1),:),'m--');
                    txttit1=['KENDALLS TAU:        ',Xfnlab,'  and  ',Yfnlab ',  ' txttau];
                    title(txttit1);
                    %pltext(.1,.9,10,[num2str(len),'-year segments']);
                    %pltext(.1,.8,10,[' Offset by ',num2str(offset),' years']);
                    ylabel('Kendalls tau');
                    xlabel('Ending year');
                    kfirst=0;
                    hh=uicontrol('Style','Pushbutton',...
                        'Units','Normalized',...
                        'Position',[.9 .9 .1 .1],...
                        'Callback','print -dps','String','Print');
                    
                    
                elseif k1==2;  % Initial and Zoom plots on ring width
                    figure(2);
                    plot(yrx,x,'-',yry,y,'--');
                    title(['RW: Solid= ',Xfnlab,'    Dashed= ',Yfnlab]);
                    xlabel('Year');
                    ylabel('RW (mm)');
                    
                    k=1;
                    while k~=4;
                        uiwait(msgbox('Click twice, first on left, then on right, to define window','Message','modal')); 
                        [tbeg,v]=ginput(1);  % get coordinates of end points of segment
                        [tend,v]=ginput(1);  % get coordinates of end points of segment
                        % Round off toward inside
                        tbeg=ceil(tbeg);
                        tend=floor(tend);
                        
                        if tbeg< max([yrx(1) yry(1)])
                            tbeg=max([yrx(1) yry(1)]);
                        end
                        if tend > min([yrx(length(x))  yry(length(y))]);
                            tend=min([yrx(length(x))  yry(length(y))]);
                        end
                        
                        tseg=(tbeg:tend)';  % time segment for zoom
                        ix = tseg-yrx(1)+1;  % index into x
                        iy = tseg-yry(1)+1;  % index into y
                        figure(2);
                        plot(tseg,x(ix),'-',tseg,y(iy),'--');
                        grid
                        title(['RW: Solid= ',Xfnlab,'    Dashed = ',Yfnlab]);
                        xlabel('Year');
                        ylabel('Ring width (mm)');
                        
                        %**********
                        
                        tseg2=tseg;
                        tseg2(1)=[];
                        iw = tseg2-yrw(1)+1;  % index into w
                        iz = tseg2-yrz(1)+1;  % index into z
                        
                        % Kendall's tau and whether signif at 95%, 99%
                        [T,tau]=kendtau(w(iw),z(iz));
                        nT=length(w(iw));  % sample size for entering lookup table
                        %Tcrit=table1(cona12,nT); % get a row from the lookup table -- OBSOLETER
                        Tcrit=interp1(cona12(:,1),cona12(:,2:size(cona12,2)),nT); % REVISED, 10-23-02
                        
                        T95=Tcrit(2);  % 95 % signif level
                        T99=Tcrit(4);  % 99% signif level
                        if T>T99,
                            sigT='**';
                        elseif T>T95,
                            sigT='*';
                        else
                            sigT=' ';
                        end
                        
                        
                        % Sign test
                        [nagree,ntot,n95,n99]=signtest(w(iw),z(iz),1); 
                        if n95=='Y' & n99=='N';
                            sig='*';
                        elseif n99=='Y';
                            sig='**';
                        else
                            sig=' ';
                        end
                        
                        
                        % Hypergeometric test, for below 0.2 quantile
                        [Nb,nsb,pb]=hypernew(w(iw),z(iz),1,.2);  % 
                        if pb<=0.05 & pb >0.01;
                            sigb='*';  % p-value 0.05 
                        elseif pb <= 0.01;
                            sigb='**';  % p-value 0.01
                        else
                            sigb=' ';
                        end 
                        
                        figure(3);
                        
                        plot(tseg2,w(iw),'-',tseg2,z(iz),'--');
                        grid;
                        title(['RW Change:  Solid= ',Xfnlab,'    Dashed = ',Yfnlab]);
                        xlabel ('Year');
                        ylabel('Change in Ring width (transformed)');
                        pltext(.5,.95,10,['tau = ',num2str(tau),' ',sigT]);
                        
                        pltext(.5,.90,10,['Sign: ',num2str(nagree),...
                                '/',num2str(ntot),' ',sig]);
                        
                        
                        pltext(.5,.85,10,['Hyperg: ',num2str(nsb),'/',...
                                num2str(Nb),' ',sigb]);
                        
                        pltext(.1,.1,10,['N = ',int2str(length(w(iw)))]);
                        figure(2);
                        
                        %***************
                        
                        
                        
                        
                        k=0;
                        while k~=1;
                            k=menu('Choose one: ','Zoom','Full','Print','Menu');
                            if k==1;
                            elseif k==2
                                plot(yrx,x,'-',yry,y,'--');
                                title(['RW: Solid= ',Xfnlab,'   Dashed= ',Yfnlab]);
                                xlabel('Year');
                                ylabel('Ring Width (mm)');
                            elseif k==3
                                eval('print');
                            elseif k==4
                                close
                                break
                            end
                        end ; % while loop on k
                    end
                    
                elseif k1==3;  %zoom analysis on transformed rw
                    %  Now transform ring width and do same
                    
                    plot(yrw,w,'-',yrz,z,'--');
                    title(['RW Change:  Solid= ',Xfnlab,'    Dashed = ',Yfnlab]);
                    xlabel('Year');
                    ylabel('Change in Ring Width (transformed)');
                    
                    k=1;
                    while k~=4;
                        uiwait(msgbox('Click twice, first on left, then on right, to define window','Message','modal')); 
                        [tbeg,v]=ginput(1);  % get coordinates of end points of segment
                        [tend,v]=ginput(1);  % get coordinates of end points of segment
                        if tbeg< max([yrw(1) yrz(1)])
                            tbeg=max([yrw(1) yrz(1)]);
                        end
                        if tend > min([yrw(length(w))  yrz(length(z))]);
                            tend=min([yrw(length(w))  yrz(length(z))]);
                        end
                        
                        tseg=(tbeg:tend)';  % time segment for zoom
                        iw = tseg-yrw(1)+1;  % index into w
                        iz = tseg-yrz(1)+1;  % index into z
                        
                        % Kendall's tau and whether signif at 95%, 99%
                        [T,tau]=kendtau(w(iw),z(iz));
                        nT=length(w(iw));  % sample size for entering lookup table
                        %Tcrit=table1(cona12,nT); % get a row from the lookup table -- OBSOLETE
                        Tcrit=interp1(cona12(:,1),cona12(:,2:size(cona12,2)),nT); % REVISED, 10-23-02
                        T95=Tcrit(2);  % 95 % signif level
                        T99=Tcrit(4);  % 99% signif level
                        if T>T99,
                            sigT='**';
                        elseif T>T95,
                            sigT='*';
                        else
                            sigT=' ';
                        end
                        
                        
                        % Sign test
                        [nagree,ntot,n95,n99]=signtest(w(iw),z(iz),1); 
                        if n95=='Y' & n99=='N';
                            sig='*';
                        elseif n99=='Y';
                            sig='**';
                        else
                            sig=' ';
                        end
                        
                        
                        % Hypergeometric test for below 0.2 quantile
                        [Nb,nsb,pb]=hypernew(w(iw),z(iz),1,.2);  
                        if pb<=0.05 & pb >0.01;
                            sigb='*';  % p-value 0.05 
                        elseif pb <= 0.01;
                            sigb='**';  % p-value 0.01
                        else
                            sigb=' ';
                        end 
                        
                        plot(tseg,w(iw),'-',tseg,z(iz),'--');
                        grid;
                        title(['RW Change:  Solid= ',Xfnlab,'    Dashed = ',Yfnlab]);
                        xlabel ('Year');
                        ylabel('Change in Ring width (transformed)');
                        pltext(.5,.95,10,['tau = ',num2str(tau),' ',sigT]);
                        
                        pltext(.5,.90,10,['Sign: ',num2str(nagree),...
                                '/',num2str(ntot),' ',sig]);
                        
                        
                        pltext(.5,.85,10,['Hyper: ',num2str(nsb),'/',...
                                num2str(Nb),' ',sigb]);
                        
                        pltext(.1,.1,10,['N = ',int2str(length(w(iw)))]);
                        
                        k=0;
                        while k~=1;
                            k=menu('Choose one: ','Zoom','Full','Print','Menu');
                            if k==1;
                            elseif k==2
                                plot(yrw,w,'-',yrz,z,'--');
                                title(['RW Change:  Solid= ',Xfnlab,'    Dashed = ',Yfnlab]);
                                xlabel('Year');
                                ylabel('RW change from previous year (scaled)')
                            elseif k==3,
                                eval('print -dps');
                            elseif k==4,
                                close
                                break
                            end
                        end;  % while loop on k
                    end 
                    
                    
                elseif k1==4; % on outer while loop (k1)
                end
            end % of while loop  (k1)
            
        elseif kk1==2; %edit ringwidths
            [X,loggy]=rwedit(X);
            E=[E; loggy];
            kfirst=1; % will want to automatically choose to analyze again
            
        elseif kk1==3; % Overwrite .rw file with edited ringwidths
            rwwrite(X,perx,whenx,Xfn)
            % also write log ascii file keeping track of editing 
            % and name aatemp.txt
            E(1,:)=[];  % first row was blanks
            fid2=fopen('aatemp.txt','wt+'); 
            fprintf(fid2,'%15s\n',E');
            fclose(fid2);
        elseif kk1==4; % make edit window
            % Get ringwidths for data listing
            tseggo = min(tseg); tsegsp=max(tseg);
            Lx1 = yrx>=tseggo & yrx<=tsegsp;
            Ly1 = yry>=tseggo & yry<=tsegsp;
            ylist=y(Ly1);
            xlist=x(Lx1);
            
            datin={tseg,x(ix),y(iy),w(iw),z(iz),...
                    yreh,vv,nH,tau95,tau99,xlist,ylist};
            figspecs={.85,.7,.05,[.1 .1],.5};
            strin={Xfnlab,Yfnlab,txttit1};
            figwind=4;
            figure(4);
            close(4);
            rwwind01(figspecs,strin,datin,figwind);
            orient tall;
            hh1=uicontrol('Style','Pushbutton',...
                'Units','Normalized',...
                'Position',[.9 .93 .1 .07],...
                'Callback','print -dps','String','Print');
            
        elseif kk1==5;
            close
            break
        elseif kk1==6;
            kk2=0;
            close;
            break;
            
        end ; % of kk1 if loop 
    end; %of kk1 while
    
end; % of kk2 while

%-- SUBFUNCTIONS

function keypress(hp1)
xdata=get(hp1,'XData');
ydata=get(hp1,'YData');
[t,x]=ginput(1);
tnear = round(t);
L=tnear==xdata;
xnear=ydata(L);
tstr=int2str(tnear);
xstr=num2str(xnear);
xlims=get(gca,'XLim');
if t>=xlims(1) & t<=xlims(2);
   text(round(t),x,{tstr,xstr},'HorizontalAlignment','center');
else;
   text(round(t),x,[]);
end;





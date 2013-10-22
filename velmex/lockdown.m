function lockdown(jpick)
% lockdown: test crossdating of dated series or find highest-correlation position of undated series
% lockdown (jpick);
% Last revised 2006-10-7
%
% Test crossdating of dated series or find highest-correlation position of undated series. For checking dated
% series, the dated series may be alternatively checked against three types of "masters": 1) a site chronology in an
% ITRDB-formatted .crn file, 2) a chronology built from all other cores in the rwl file, or 3) a chronology
% built from a subset of ring-width series whose names are read from an input .mat file prepared previously.
%
%
%*** INPUT
%
% jpick (1 x ?)i , [], or missing:  pointer to series in the rwl file to be tested in this run;  allow for the analysis to
%     be restricted to a few test cores, which might be necessary with some huge datasets (e.g., hundreds of cores).
%     Example 1:  jpick=[1 3 5] use only the first, third and fifth core in the rwl file
%     Example 2:  jpick=[34:end]   use only the 34th through last cores
%     Example 3:  jpick =[] use all cores.  Missing argument jpick also means this
%.
% User prompted for the name of the .rwl file with the ring widths to be tested, and optionally prompted for
% names of a .crn file with a master chronology or a .mat file with core ids of series for the tentative
% master.
%
%*** OUTPUT
%
% Screen plots summarizing results:
% 1) line plot with minimum, maximum and median r of each ring-width series and the master
% 2) Grayscale map flagging for each ring-width series and m1-year period offset by m2 years whether
%   1) the highest r is at the dated position and r>rc, where rc is a specified critical threshold correlation
%   2) the highest r is at the dated position, but r<=rc
%   3) the highest r is not at the dated position, but r>rc at the dated position
%   4) the highest r is not at the dated position, and r<=rc at the dated position
% 3) Colormap similar to map in (2), but going from red (highest r) to blue (lowest r)
% 4) Time series plot windowed r on time (top), standard devs (middle), and with the corresponding
%   high-pass z-score series in zoomable window (bottom)
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED
%
%*** TOOLBOXES NEEDED
%
%*** NOTES
%
%
% Correlations are computed with a window width m1 years and offset of m2 years.  The series are organized
% so that the first correlatons is for times 1:m1, second is for (1+m2):(m1+m2), etc.  If the length of the
% series is not be divisible by m, a final m-year period is aligned to end with last year of the common period.
% Default changeable m1 and m2 are hardcoded. Reasonable choices are m1=50 and m2=5;
%
% Series whose dating is in question are called "test" series.  The series used to check the dating of the test
% series is called the "reference" series
%
% To emphasize high-frequency ring-width variation, which must be relied on for accurate annual dating, series are
% high-pass filtered with a cubic smoothing spline g and then converted to zscores.  The high-pass version of a series
% is computed by the ratio method ( v = w/g) if dealing with core ring-widths from the rwl file, and by the
% difference method (v = w-g) if w is a site chronology from a crn file.  The master is computed by arithmetic
% averaging over series after they have been individually converted to zscores.  The average series itself is then
% converted to a zscore.  Thus the master, by whatever method, has mean 0 and std 1.
%
%
% If the master series is built from ring widths of a subset of cores, the following operations are taken:
% 1) each ring-width series is high-pass-filtered and then converted to unit variance, as described above
% 2) the transformed ring-width series are averaged together to form a mean time series
% 3) that mean series is again converted to unit variance
%
% The cubic smoothing spline is specified in an input dialog by the wavelength at which the frequency response
% is 0.5.  A reasonable choice for this wavelength is 20 years
%
% The threshold critical correlation is set to r_c = 2/sqrt(m1), where m1 is the window width
%
% Test and master subsets of ringwidths from the same rwl file are kept mutually exclusive.  In other words, no series is ever
% tested against a master that includes that series.  If the master is a crn file, there may be overlap in test
% and master series.  For example, the crn may include cores also being tested.
%
% The master is not changed during a run of lockdown.  A single master applies, except if kref==2 a set of n
% master series applies, where n is the number of test cores and for each test core the master is an average
% over all other cores. 
%
% In the standard mode, sample size for the master is annotated above the x axis at years evenly divisible by
% 20 and at the first year of overlap of master and test series.  To avoid cluttered text, the first annotated 20-year sample size is 
% suppressed if it falls less then 10 years after the first year of overlap of test and master. 
%
%
% If master is built using option for "all other series", the all other is limited to cores specified by 
% input argument jpick.  For example, if jpick =[] or the jpick argument is not include, "all other" means all
% other cores in the rwl file containing the test cores.  But, say jpick==1:10.  Then core series 1,2,...10 are
% the only cores tested, and the only cores considered for the master.  When testing series 1, the master is
% made up of series 2 -10.  When testing series 2, the master is made from the set {1 3-10}, etc.
%
% Notation
% variable 1: test series
% variable 2: reference series
% u full length original series
% v windowed original series
% w windowed scaled high-pass series
% s standard-deviation of rw change, test series

close all;
clc;

if nargin==0;
    jpick=[];
end


%--- HARD CODE

wind_width=50; % default window width
rthresh = 2/sqrt(wind_width); % critical threshold
wind_shift=5; % initialize window width
spline_period=20; % wavelength spline has frequency response 0.
p_spline = splinep(spline_period,0.50); % spline parameter
str_spline1 = [num2str(spline_period) ' spline'];
str_spline2 = ['High Pass (' num2str(spline_period) ' spline)'];



%--- CHOOSE TO DEAL WITH EITHER TENTATIVELY DATED SERIES OR 'UNDATED' SERIES DEFINED TO HAVE STARTING YEAR 1
kmode = menu('Choose mode','Test tentatively dated series','Blind shifting of a test series against master')';


if kmode==1;


    %--- CHOOSE FORM OF REFERENCE SERIES: CRN FILE, ALL OTHER CORE RING WIDTHS, SPECIFED SUBSET OF CORE RING WIDTHS

    kref = menu('Choose form of master','A crn-file site chronology',...
        'Scaled high-pass average of all other core indices',...
        'Scaled high-pass average of core indices whose ids will be read in');



    %--- OPTIONALLY READ A MASTER CRN FILE

    if kref==1;
        [file2,path2]=uigetfile('*.crn','crn Infile with master series');
        pf2=[path2 file2]; % merge path and filename
        [y,s,yry]=crn2vec2(pf2); % y is a site chronology, s is sample size, yry is year
        y1 =y;
        yry1=yry;
        s1=s;
        [y1,yry1]=trimnan(y1,yry1);
        [s1,yrs1]=trimnan(s1,yry);
    else;
        % No action; will be building master from the rwl file
    end;


    %--- INPUT DIALOG TO SPECIFY SPINE WAVELENGTH, WINDOW WIDTH

    prompt={'Enter the spline wavelength:','Enter the sliding window width:','Enter the offset:'};
    name='Program settings';
    numlines=1;
    defaultanswer={num2str(spline_period),num2str(wind_width),num2str(wind_shift)};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    spline_period=str2num(answer{1});
    wind_width=str2num(answer{2});
    wind_shift = str2num(answer{3});
    p_spline = splinep(spline_period,0.50); % spline parameter



    %--- INPUT THE RWL FILE (HAS THE TEST SERIES, AND MAYBE THE SERIES FOR THE MASTER)

    [file1,path1]=uigetfile('*.rwl','rwl infile with ''test'' series');
    pf1=[path1 file1]; % merge path and filename
    pf1a=rwlinp(pf1);
    eval(['load ' pf1a ';']); % X, nms, yrs are key
    if size(jpick,1)==1;
        jpick=jpick';
    end;
    if isempty(jpick);
        jpick=(1:size(nms,1))';
    end
    Xall=X;
    nmsall=nms;
    yrsall=yrs;
    jall = (1:size(nmsall,1))';


    %--- OPTIONALLY READ A .MAT FILE OF CORE IDS OF MASTER SERIES

    if kref ==3;
        [file_ids,path_ids]=uigetfile('ids_*.mat','Infile with ids of cores for building master');
        pf_ids =[path_ids file_ids];
        eval(['load ' pf_ids ';']);
        if ~exist('ids_master','var');
            error([pf_ids ' does not have variable ids_master']);
        else;
            if ~isa(ids_master,'cell');
                error('ids_master not a cell');
            end
        end
        nref = length(ids_master); % number of cores in specified master
        file2=file1;
        
        % Update jpick by omitting any master cores -- do not want to test any of the master cores
        dtemp=cellstr(nmsall);
        jmaster = find(ismember(nmsall,ids_master));
        jpick = setxor(jpick,jmaster);
        if isempty(jpick);
            error('set of test cores is empty after masking those in master');
        end

        % Build a string menu of ids of cores, with those in master marked M
        nmsall_cell=cellstr(nmsall);
        nmstest_cell=nmsall_cell(jpick);
        nmsmaster_cell =nmsall_cell(jmaster);
        L=ismember(nmsall_cell,nmsmaster_cell) | ismember(nmsall_cell,nmstest_cell);
    elseif kref==1; % using crn file as reference
        nmsall_cell=cellstr(nmsall);
        nmstest_cell=nmsall_cell(jpick);
        nmsmaster_cell=[];
        pf_ids = [path1 'ids_master.mat'];
    elseif kref==2; % using all other cores as master
        nmsall_cell=cellstr(nmsall);
        nmstest_cell=nmsall_cell(jpick);
        file2=file1;
        pf_ids = [path1 'ids_master_initial.mat'];
    end

    kmaster=0; % master has not been changed yet


    %--- STORE TEST SERIES IN TIME SERIES MATRIX

    [X,yrX,nms]=sov2tsm3(X,yrs,nms,jpick,[]); % x is tsm of rw series, yrx is year vector; nms is col-cell of strings
    [mX,nX]=size(X);
    ntest = nX; % this many ring width series to be tested

    % Check that the years of selected test series are consistent with the dated-series  mode
    T = firstlst(X,yrX);
    if any(T(:,1)==1);
        error(['You said dated-series mode, but one or more series starts with year 1']);
    end



    %--- COMPUTE ZSCORE HIGH PASS MASTER
    %  If just one master for the analysis (kref==1 or kref==3), will have it in
    % x_master, yrx_master,n_master, where n_master is the time-varying sample size for master
    % If more than one master (kref==2), have X_master, yrX_master, SS,yrSS, where these are time series
    % matrices the same col-size as the number of test series, and each column is the master to be used for
    % that test series.

    % Store all ring widths in time series matrix
    [Xall,yrXall,nmsall]=sov2tsm3(Xall,yrsall,nmsall,((1:size(nmsall,1))'),[]);
    [mXall,nXall]=size(Xall); % the "all" variables may be later needed if updating the master

    if kref==1; % using crn file for master
        % Recall have crn master in y1, yry1, with sample size in s1
        n_master = s1;
        kdir = 2; % use difference to compute high pass
        [x_master,yrx_master]=subfun01(y1,yry1,p_spline,kdir);
        
        % Build time series of master sample size (number of cores or trees, depending on source crn)
        yrss =yry1;
        ss = s; 
             
    elseif kref==2; % using all non-test cores for master -- will have a different master for each test core

        %--- Build master
        %
        % The master, for this option, is the average s-score high-pass index for "all other" cores -- the set
        % not included the core tested.  It is assumed that all cores no        t in jpick are NOT to be considered to
        % be tested or to be used in building a master
        
        M = repmat(NaN,mX,nX); % to hold the master series (zscore high-pass) applicable to each core
        yrM=yrX;
        A = repmat(NaN,mX,nX); % to hold zscore high-pass series (for each core)
        yrA=yrX;
        G = repmat(NaN,mX,nX); % to hold spline curve used to high pass each series
        yrG=yrX;

        % Check that the years of selected test series are consistent with the dated-series  mode
        T = firstlst(X,yrX);
        if any(T(:,1)==1);
            error(['You said dated-series mode, but one or more series starts with year 1']);
        end
        YRS1 = T; % store first and last valid year of each test series
        
        % Generate zscore high-pass versions of all jpick series
        for n =1:nX; % loop over all seris in jpick
            x = X(:,n);
            yrx = yrX;
            kdir=1; % use ratio in computing high-pass
            [y,yry,gv]=subfun01(x,yrx,p_spline,kdir);
            igo = yry(1)-yrM(1)+1;
            isp =igo + length(y) -1;
            A(igo:isp,n)=y; % store z-score high pass test series
            G(igo:isp,n)=gv; % store spline curves used in computing A
        end; %  for n =1:nX; % loop over all seris in jpick

        
        
        %----  Compute the nX master series
        
        nmaster=nX-1; % number of series available (not necessarily all in same year ) for each master
        M = repmat(NaN,mX,nX);  % to hold each test-specific master sieres
        LLmaster = logical(zeros(nXall,ntest));  % to hold pointer to cores in each master in Xall, nmsALL
        Nms_master = cell(nmaster,ntest); % each col of this cell-mtx has the ids of cores in the  a master
        SS = zeros(mX,nX); % to hold time-varying sample depths of masters
        yrSS=yrX;
        
        
      
        Lshort = logical(zeros(nX,1)); % to flag test series with insufficient overlap with reference series
        L = logical(ones(nX,1));
        for n =1:nX; % loop over columns of A
            nmsthis = nms;
            nmsthis(n)=[];
            Nms_master(:,n) = nmsthis;
            L1 = L;
            L1(n)=0;
            A1 = A(:,L1);
            yrA1=yrA;
            
            T = firstlst(A1,yrA1); % 2-col mtx of first and last years of series
            ton = min(T(:,1));
            toff = max(T(:,2));
            L4= yrA>=ton & yrA<=toff;
            igo = ton - yrA(1)+1;
            isp = toff-yrA(1)+1;
            A1 = A1(L4,:);
            yrA1=yrA1(L4);
            [mtemp,ntemp]=size(A1);
            if ntemp==1;
                ymaster=A1;
            else
                ymaster = (nanmean(A1'))'; % average over cores to be used as master
            end
            M(igo:isp,n)= ymaster; % Store masters

            % Convert master to zscores
            mn_temp=nanmean(ymaster);
            sd_temp=nanstd(ymaster);
            ymaster = (y - mn_temp) / sd_temp;
            clear mn_temp sd_temp;
            
            % Compute time-varying sample size (number of cores) in masters
            Z = (A1)';
            L5 = ~isnan(Z);
            n_master = (sum(L5))';
            SS(igo:isp,n)=n_master;
            clear n_master L5 Z A1 yrA1 toff ton igo isp T
        end; % for n =1:nX; % loop over columns of A
        YRS2 = firstlst(M,yrM); % first and last year of each master series
        
        
        % Flag series with insufficient overlap with master for even 1 windowed  correlation
        for n = 1:ntest; % loop over test series
            yruz = (YRS1(n,1):YRS1(n,2))'; % year cv for test series
            yrvz = (YRS2(n,1):YRS2(n,2))'; % year cv for master
            nlap = length(intersect(yruz,yrvz));
            if nlap<wind_width;
                Lshort(n)=1;
            end
        end

              
        clear L yruz yrvz n nlap
        % Now have this:
        %
        % M, yrM;  master series, a time series matrix, same size as X
        % YRS1: 2-col matrix of first and last valid year of each master
        % SS, yrSS:   sample size for M
        %  Nms_master: a cell-string matrix, each col with the ids of cores in the master for that col of X
        % Lshort: cv with zeros, except 1 if less than wind_width overlap between master and test

    elseif kref==3; % building master from selected cores
        Lmaster = ismember(nmsall,ids_master);
        nmaster = length(ids_master);
        if sum(Lmaster) ~= nmaster;
            error(['You specified ' num2str(nmaster) ' cores as masters, but only ' sum(Lmaster) ' name matches of ids']);
        end
        X_master = Xall(:,Lmaster);
        yrX_master = yrXall;
        nms_master=ids_master;
        Y_master = repmat(NaN,size(X_master,1),nmaster); % to hold zscore high-pass version of X_master
        Gv = Y_master; % for spline curves used in coverting X_master to Y_master
        yrY_master = yrX_master;
        for n = 1:nmaster; % loop over master series, computing zscore of high-pass
            x = X_master(:,n);
            yrx = yrX_master;
            kdir=1; % use ratio in computing high-pass
            [y,yry,gv]=subfun01(x,yrx,p_spline,kdir);
            igo = yry(1)-yrY_master(1)+1;
            isp =igo + length(y) -1;
            Y_master(igo:isp,n)=y;
            Gv(igo:isp,n)=gv;
        end;

        % Average over master cores
        [mtemp,ntemp]=size(Y_master);
        if ntemp>1;
            T = firstlst(Y_master,yrY_master); % 2-col mtx of first and last years of series
            ton = min(T(:,1));
            toff = max(T(:,2));
            L= yrY_master>=ton & yrY_master<=toff;
            Y_master = Y_master(L,:);
            yrY_master=yrY_master(L);
            Z = (Y_master)';
            L = ~isnan(Z);
            n_master = (sum(L))';
            
            % Make time series of sample size (number of cores) in master
            yrss = yrY_master; 
            ss = n_master;
            
            z = (nanmean(Z))'; % column vector of mean of the zscores of individual high-pass master cores
            yrx_master = yrY_master;
            x_master = (z - nanmean(z)) / nanstd(z); % re-scale average series as zscore
        else
            x_master=y;
            yrx_master = yry;
            n_master = ones(length(y),1); % number of cores in master
            
            % Make time series of sample size (number of cores) in master
            yrss = yrx_master; 
            ss = n_master;
        end
    else; % all-other cores master
        error('Must code for sample size time series of master for all-other-series option');
    end
    
    
    %--- HONE THE PLOTTING TIME SERIES OF SAMPLE SIZE FOR MASTER
    %
    % On plot, will want to annotate sample size at first year and at every year evenly divisible
    %  by 20
    
    if kref~=2;
        L1 = mod(yrss,20)==0;
        L1(1)=1;
        yrss1=yrss(L1);
        ss1=ss(L1);
    else;  % "all others" are masters;

        % Make cells with time series of sample size for masters
        sscell = cell(ntest,1);
        for n=1:ntest;
            ss = SS(:,n);
            yrss=yrSS;
            L2 = ss>0;
            ss=ss(L2);
            yrss=yrss(L2);
            if ~all(diff(yrss)==1);
                error('yrss not inc by 1');
            end
            L1 = mod(yrss,20)==0;
            yrss1=yrss(L1);
            ss1=ss(L1);
            sscell{n}=[yrss1 ss1];
        end
    end


    %--- SPLINE-FIT, HIGH-PASS FILTER AND SCALE ALL TEST  SERIES; ALSO CHECK OVERLAP WITH MASTER
    %
    % Recall tha X, yrX, nms has the time series, year and names, and that this set has already been culled
    % according to jpick.
    
    
    if kref~=2; % Recall that if kref==2, have already done the high-pass filtering and computed the master
        %   series (multiple in that case), as well as other matrices, such as A, YRS1, G, Lshort, etc

        A = repmat(NaN,mX,nX); % to store the scaled high-pass test series
        YRS1 = repmat(NaN,nX,2); % to store first and last year of each test series
        G = A; % to store growth curves
        yrA= yrX;
        Lshort = logical(zeros(nX,1)); % to flag test series with insufficient overlap with reference series
        for n = 1:nX; % loop over test series
            if kref ==1 | kref==3; % have a single relevant master
                vz =x_master;
                yrvz =yrx_master;
            else;
                vz = X_master(:,n);
                yrvz = yrX_master;
                [vz,yrvz]=trimnan(vz,yrvz);
            end
            x = X(:,n);
            yrx = yrX;
            kdir=1; % use ratio in computing high-pass
            [uz,yruz,gu]=subfun01(x,yrx,p_spline,kdir);
            YRS1(n,:)=[yruz(1) yruz(end)];
            igo = yruz(1)-yrA(1)+1;
            isp =igo + length(uz) -1;
            A(igo:isp,n)=uz;
            G(igo:isp,n)=gu;
            % Check for overlap with vz, yrvz
            nlap = length(intersect(yruz,yrvz));
            if nlap<(wind_width+2); % extra 2 yrs overlap to allow for computing lagged r
                Lshort(n)=1;
            end
        end; % for n = 1:nX; % loop over test series
        clear n vz yrvz x yrx uz yruz gu igo isp nlap;
    else
    end
    
    % Dead in water if no series have enough overlap with master
    if sum(~Lshort)==0;
        error(['No test series marked by jpick have at least ' num2str(wind_width) ' yr overlap with reference']);
    end
    
    
    %---- Warn if some test series too short for testing with master
    nshort = sum(Lshort); % number of test series with insufficient overlap with master
    if nshort>0;
        nmsout =nmstest_cell(Lshort);
        mess1 = ['These series have < ' num2str(wind_width+2) ' years overlap with master--will not be tested'];
        mess1=[mess1; nmsout];
        uiwait(msgbox(mess1,'Message','modal'));
    end
   
    
    %--- BUILD CHAR MATRIX WITH FIRST AND LAST YEAR OF ALL CORE SERIES FROM THE RWL (TEST, MASTER, AND ANY NOT
    % IN JPICK

    str_men2=[repmat('(',nXall,1) num2str(yrsall(:,1)) repmat('-',nXall,1)  num2str(yrsall(:,2)) repmat(')',nXall,1)];
    menu_test2 = [num2str(jall) repmat(' ',nXall,1) char(nmsall) repmat(' ',nXall,1) str_men2];
    menu2 = cellstr(menu_test2);
    clear str_men2 menu_test2;
    
    
    %--- ELIMINATE THE TEST SERIES WITH TOO-SHORT OVERLAP WITH MASTER
    
    jpick(Lshort)=[];
    X = X(:,~Lshort);
    nX = size(X,2);
    G(:,Lshort)=[];
    A(:,Lshort)=[];
    YRS1 = YRS1(~Lshort,:);
    nms = nms(~Lshort);
    ntest=length(nms); 
    
    % Depending, also eliminate the masters
    if kref==2 & any(Lshort);
        M(:,Lshort)=[];
        [mM,nM]=size(M);
        Nms_master(:,Lshort)=[];
        sscell(Lshort)=[];
    end;
    

    %--- BUILD MENU CHOICES TO REMAINING TEST SERIES
     
    str_men1=[repmat('(',nX,1) num2str(YRS1(:,1)) repmat('-',nX,1)  num2str(YRS1(:,2)) repmat(')',nX,1)];
    menu_test = [num2str(jpick) repmat(' ',nX,1) char(nms) repmat(' ',nX,1) str_men1];
    menu1 = cellstr(menu_test);
    clear str_men1 menu_test;


    %---- COMPUTE AND STORE THE CORRELATIONS AND INDIVIDUAL-PLOT INFORMATION

    % The three subplots with time series plots of windowed r, standard deviations, and z-score series
    F1_titles = cell(ntest,1); % title
    F1_legends = cell(ntest,1); % legends
    F1_r = cell(ntest,1); % 4-col tsm of end year , correlations of test_t vs master_t, test_t vs master_t+1,
    %      test_t vs master_t-1
    F1_s = cell(ntest,1); % 3-col tsm of year, windowed std of test, of master


    % The summary correlation figure
    F2_legends = {'Median','Minimum','Maximum'};
    F2_xlabel='Series Number';
    F2_ylabel='Correlation';
    F2_xticklabel=cellstr(num2str(jpick));
    F2_data = repmat(NaN,ntest,3);  % to hold median, minimum and maximum correlation

    str_lag = {'x_{t} vs y_{t}','x_{t} vs y_{t+1}','x_{t} vs y_{t-1}'}; % this will be a legend on top plot

    for n = 1:ntest; % loop over test series
        nm_ref = strtok(file2,'.'); % name of master, or reference, series
        nm_test = [strtok(file1,'.') '-' nms{n}]; % name of test series
        if kref==1;
            nm_ref = strtok(file2,'.'); % name of ref series
        elseif kref==2;
            nm_ref = [strtok(file1,'.') '- all others'];
        else
            nm_ref = [strtok(file1,'.') '- specified'];
        end
        F1_legends{n}={nm_ref,nms{n}};

        % Get pair of series
        x = A(:,n); % test series, z score high pass
        yrx=yrA;
        if kref==1 | kref==3;
            y = x_master;
            yry = yrx_master;
        else; % master is made from all other series in jpick 
            y=M(:,n);
            yry = yrM;
        end

        % Trim
        [x,yrx]=trimnan(x,yrx);
        [y,yry]=trimnan(y,yry);
        str_full=[num2str(yrx(1)) '-' num2str(yrx(end))]; % full coverage of test series


        % FIND OVERLAP OF HIGH-PASS AND CULL

        yrgo = max([yrx(1) yry(1)]);
        yrsp = min([yrx(end) yry(end)]);
        str_checkable=[num2str(yrgo) '-' num2str(yrsp)]; % full coverage of test series
        L=yrx>=yrgo & yrx<=yrsp;
        x=x(L);
        yrx = yrx(L);
        L=yry>=yrgo & yry<=yrsp;
        y=y(L);
        yry = yry(L);

        if sum(L)<(wind_width+2);
            str_hey1 = [nms{n} ' (Series # ' num2str(jpick(n)) ' in ' file1 ') short; maximum allowable window width is ' num2str(sum(L)-1) ' yr'];
            error(str_hey1);
        end

        str_1=['All: ' str_full '; Checkable: ' str_checkable];
        F1_titles{n}=['Test series, x_{t} = ' nms{n} ' (' str_1 ');   Master series, y_{t} = ' nm_ref];

        % STORE CULLED HIGH-PASS DATA FOR LATER INTERACTIVE PLOTS
        D.data{n}=[yrx x y];
        D.nmtest{n}  = nms{n};
        D.nmref{n} = nm_ref;

        % COMPUTE DIFFERENCE OF HIGH-PASS TEST AND HIGH-PASS REF
        % Want to compute absolute departures and get their 0.95 prob point.  Then build a logical pointer to years
        % in which absolute departure is unusually large.  Later, will compute sums of those large-values in the
        % moving window, and plot at end of window year.
        d = x-y;
        dabs=abs(d);
        Lbig=dabs>prctile(dabs,95);


        % MAKE MATRICES OF INDEX POINTERS TO STARTING AND ENDING YEARS OF WINDOW SEGMENTS--I3

        nsize = length(x);
        yrend = ((yrx(1)+wind_width-1):wind_shift:yrx(end))';
        if yrend(end) < yrx(end);
            yrend=[yrend; yrx(end)];
            yrstart = yrend-wind_width+1;
        end;
        D.yrend{n}=yrend;
        isp = (yrend-yrx(1)+1)'; % rv of ending rows of x and y for each segment
        I1 = repmat(isp,wind_width,1);
        i2 = flipud((0:(wind_width-1))');
        I2 = repmat(i2,1,size(I1,2));
        I3 = I1-I2;
        clear isp I1 I2 i2;

        % Shifted pointers; X(I3fa) vs Y(I3ra) is test in year t with ref in t+1
        %  X(I3fb) vs Y(I3rb) is test in year t with ref in t-1
        I3fa = I3(1:(end-1),:);
        I3ra = I3((2:end),:);
        I3fb = I3((2:end),:);
        I3rb = I3((1:end-1),:);

        
        % MAKE TSMs OF SHIFTED TIME SERIES OF  SERIES -- F and R for test and reference

        F = x(I3);
        R = y(I3);
        B = Lbig(I3);


        %-- STANDARD DEVIATIONS, WINDOWED

        s_f = nanstd(F);
        s_r= nanstd(R);
        F1_s{n}=[yrend s_f' s_r'];


        %--- CALL CORRPAIR TO GET THE CORRELATIONS FOR EACH SUBPERIOD

        [r,nr]=corrpair(F,R);
        D.r{n}=r;
        D.summary(n,:)=[length(r) median(r) min(r) max(r)];
        nbig = nansum(B); % number of unusually large departures in zscore high-pass series in the window
        [ra,nra]=corrpair(x(I3fa),y(I3ra));
        [rb,nrb]=corrpair(x(I3fb),y(I3rb));

        %-- STORE
        F1_r{n} = [yrend r' ra' rb']; % year, correlation simultaneous, lagged -1, +1


    end; % for n = 1:ntest; % loop over test series



    %---- PLOTTING
    
    kvirgin = 1; % logical tells whether initial master has been modified yet. When kvirgin==0, it has.

    kwh1 = 1;
    while kwh1==1;

        kmen1 = menu('Choose','Check a ''test'' series','View summary stats',...
            'Modify makeup of master','Quit');
        if kmen1==1;
            kmen2 = menu('Choose test series',menu1);
            knew=1;  % you have selected another test series
            
            % May need to re-set the time-varying sample size for master
            if kref==2;
                ss1=sscell{kmen2}(:,2); 
                yrss1=sscell{kmen2}(:,1); 
                
            else
            end
                       
            kwh4=1;
            while kwh4==1
                figure(2);
                [cL,cB,cW,cH]=figsize(0.8,0.6);
                set(gcf,'Position',[cL cB cW cH]);
                if knew==1;
                    yrx = D.data{kmen2}(:,1); % the zscore high pass data, year
                    yry=yrx;
                    x=D.data{kmen2}(:,2); % test series
                    y=D.data{kmen2}(:,3); % master series

                    xlims= [yrx(1)-3 yrx(end)+3];
                    Fthis=F1_r{kmen2};
                    yrend = Fthis(:,1);
                    r = Fthis(:,2);
                    ra = Fthis(:,3);
                    rb = Fthis(:,4);
                    Fs = F1_s{kmen2};
                    leg_s = F1_legends{kmen2};

                    subplot(3,1,1);
                    hp1=plot(yrend,r,yrend,ra,yrend,rb);
                    set(hp1(1),'LineWidth',2);
                    set(gca,'XLim',xlims);
                    xlabel(['End Year of ' num2str(wind_width) '-yr window']);
                    ylabel('r');
                    hleg1 = legend(str_lag);
                    set(hleg1,'Position',[0.151    0.7850    0.0845    0.1280])
                    grid on;
                    title(F1_titles{kmen2});
                    htext1=text(yrend(1),r(1),'o','FontSize',14,'Visible','Off'); % Invisible circle to be moved later
                    
                    
                    % Subplot with standard deviations of test and master in sliding window; also the sample
                    % size (number of cores) annotated above x axis.  The sample size may be number of trees if
                    % the master is a crn file and the sample size is defined as number of trees in that file
                    subplot(3,1,2);
                    hp2=plot(Fs(:,1),Fs(:,2),Fs(:,1),Fs(:,3));
                    legend(leg_s);
                    set(gca,'XLim',xlims);
                    ylims_ss = get(gca,'YLim');
                    ypt_ss = ylims_ss(1) + 0.02*diff(ylims_ss);
                    Ltemp = yrss1>=yrx(1) & yrss1<=yrx(end);
                    ss2 = ss1(Ltemp);
                    yrss2 = yrss1(Ltemp);
                    if yrss2(1)-yrx(1)<10;
                        yrss2(1)=[];
                        ss2(1)=[];
                    else
                    end
                    yrss2 = [yrx(1); yrss2];
                    ss2 = [ ss(yrss==yrx(1)); ss2];
                    txt_ss = cellstr(num2str(ss2));
                    ypt_ss = repmat(ypt_ss,length(yrss2),1);
                    text_ss=text(yrss2,ypt_ss,txt_ss,'Horizontalalignment','Center','VerticalAlignment','Bottom');
                    xlabel(['End Year of ' num2str(wind_width) '-yr window']);
                    ylabel('Standard Dev');
                    grid on;

                    subplot(3,1,3);
                    hp3=plot(yry,y,yrx,x);
                    set(gca,'XLim',xlims);
                    D.xlims{n}=xlims;
                    D.ylims{n} = get(gca,'YLim');
                    ylabel('z-score');
                    xlabel('Year');
                    legend(leg_s);
                    grid on;
                    knew=0;
                else
                end

                kmen5 = menu('Choose','Click on a point in top plot for zoomed time series in bottom','No zoom, move on');
                kfirst = 1; % first time through, no need to replot time series yet

                if kmen5==1;

                    if kfirst==1;
                    else;
                        figure(2)
                        subplot(3,1,3);
                        hp3=plot(yry,y,yrx,x);
                        set(gca,'XLim',xlims);
                        D.xlims{n}=xlims;
                        D.ylims{n} = get(gca,'YLim');
                        ylabel('z-score');
                        xlabel('Year');
                        legend(leg_s);
                        grid on;

                    end


                    % Get end point for zoom from top series
                    figure(2);
                    subplot(3,1,1);
                    hkids=get(gca,'Children');
                    set(hkids(1),'Visible','off');
                    [xpt,ypt]=ginput(1);
                    diffx = abs(xpt-yrend);
                    [w,iw]=min(diffx);
                    yrlast = yrend(iw);
                    %htext_circ= text(yrlast,rthis(iw),'o','HorizontalAlignment','Center','VerticalAlignment','Middle');
                    set(hkids(1),'Visible','On','Position',[yrlast r(iw) 0],'HorizontalAlignment','Center','VerticalAlignment','Middle');
                    yrfirst = yrlast-wind_width+1;
                    yrslab=(yrfirst:yrlast)';
                    Lslab = yrx>=yrfirst & yrx<=yrlast;

                    subplot(3,1,3);
                    plot(yrslab,y(Lslab),'-o',yrslab,x(Lslab),'-^');
                    ylabel('z-score');
                    xlabel('Year');
                    legend(leg_s);
                    grid on;
                    kfirst =0;
                else
                    kwh4=0;
                end

            end; % while kwh4==1
        elseif kmen1==2;
            figure(1);
            ylims_temp = [min(D.summary(:,3))-0.1  1];
            if ylims_temp(1)>0;
                ylims_temp(1)=0;
            end
            jseq = (1:nX)';
            hp1 = plot(jseq,D.summary(:,2),'o','Linewidth',2,'MarkerSize',12);
            hold on;
            hp2=plot(jseq,D.summary(:,3),'x','Color',[1 0 0]);
            hp3=plot(jseq,D.summary(:,4),'x','Color',[0 0.6 0]);
            set(hp2,'MarkerSize',10,'LineWidth',2);
            set(hp3,'MarkerSize',10,'LineWidth',2);

            xlabel('Series Number');

            ylabel('Correlation');
            legend('Median','Minimum','Maximum');
            set(gca,'XTick',jseq,'XLim',[0 nX+1],'YLim',ylims_temp,'Xticklabel',num2str(jpick));
            grid on;
            hold off;
            fheight=0.3;
            fwidth=0.1 + 0.7*(nX/50);
            fwidth=min([fwidth 0.8]);
            [cL,cB,cW,cH]=figsize(fwidth,fheight);
            set(gcf,'Position',[cL cB cW cH]);
        elseif kmen1==3 ;  % Modify makeup of master
            if kref==1 | kref==3;
                if isempty(nmsmaster_cell) &   kref==1;
                    Lin = logical(zeros(length(nmsall),1)); % none yet in master
                else
                    Lin = ismember(nmsall,nmsmaster_cell);
                end
                nallow =[0 length(Lin)];
                strmenu='Toggle Y for a core to be in master, N for a core to not be in master';
                C1=menu2;
                
                C2=cellstr(repmat('-N',length(C1),1));
                C2(Lin)=cellstr('-Y');
                [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu);
                nmsmaster_cell=nmsall(Lout);
                %nmsmaster_cell=C3';
                kmaster=1;
            else; % assumes "all other cores"  in master;  here you can build a "separate" master one by
                % one by adding cores 
                C1=menu2;
                if kvirgin==1;
                    C2=cellstr(repmat('-N',length(C1),1));
                    Lin = logical(zeros(length(nmsall),1)); % none yet in master
                    kvirgin=0; % no longer virgin
                else
                    Lin = ismember(nmsall,nmsmaster_cell);
                    C2(Lin)=cellstr('-Y'); % toggle 
                end;
%                 nmsmaster_cell=Nms_master(:,kmen2);
%                 if isempty(nmsmaster_cell);
%                     Lin = logical(zeros(length(nmsall),1)); % none yet in master
%                 else
%                     Lin = ismember(nmsall,nmsmaster_cell);
%                 end
               
                strmenu='Toggle Y for a core to be in master, N for a core to not be in master';
                
                nallow =[0 length(Lin)];
                [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu);
                nmsmaster_cell=nmsall(Lout);
                %nmsmaster_cell=C3';
                kmaster=1;
                
            end

        else; % kmen1 == 4
            kwh1=0;
            

            % Optionally update the master
            if kmaster==1;
                ktemp = questdlg(['Update ' pf_ids ' with revised list of master cores?']);
                if strcmp(ktemp,'Yes');
                    ids_master = nmsmaster_cell;
                    eval(['save ' pf_ids ' ids_master;']);
                end
            else
            end

            kclose = questdlg('Close all windows?');
            if strcmp(kclose,'Yes');
                close all hidden;
            else
            end
        end;

    end; % while kwh1==1;
else; %  kmode==2
    shiftcor(jpick); % function that deals with possible dating position of undated
end; % if kmode==1;


% return  here if change the composition of master
%--- GENERATE THE HIGH-PASS MASTER

%--- COMPUTE THE INDIVIDUAL WINDOWED CORRELATIONS AND THEIR TIME SERIES PLOTS (length(jpick) windows)

%--- COMPUTER THE MATRIX OF CORRELATIONS FOR THE GRAY-SCALE AND COLOR MAPS

%--- MAKE THE GRAYSCALE AND COLORMAP FIGURES

%---- MAKE THE SUMMARY CORRELATION FIGURE

%---- MAKE THE MENU FOR
%   Check a test series
%   View summary correlation plot
%   View grayscale map of segment correlations
%   View color map of segment correlations
%   (optional--if master not crn) Change composition of master


%---- SUBFUCTIONS

function [y,yry,g]=subfun01(x,yrx,p_spline,kdir);
%
% Compute z-score high-pass of a time series
%
xtemp=x;
[xtemp,yrxtemp]=trimnan(x,yrx); % trrim leading and trailing NaNs
L=isnan(xtemp); % mark any internal NaN
xxtemp=xtemp(~L); % pull non-NaNs of trimmed
yrxxtemp=yrxtemp(~L);
[g,p] = csaps(yrxxtemp,xxtemp,p_spline,yrxtemp); % compute smoothing spline
if kdir==1; % if high pass if ratio of time series to smooth curve
    if any(g<=0);
        error('ratio high pass blows up because g goes to zero or negative');
    end
    u = xtemp ./ g; % high pass series, as ratio
else;
    u = xtemp - g; % high pass series, as difference

end
umean = nanmean(u);
ustd = nanstd(u);
y = (u-umean) / ustd; % zscore of high-pass
yry= yrxtemp;



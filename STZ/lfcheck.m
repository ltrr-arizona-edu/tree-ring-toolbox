function twarn=lfcheck
% lfcheck:  check file of ring-width series for consistency of variations as function of wavelength and detrending choice
% lfcheck; 
% Last revised 2-20-02
%
% Check file of ring-width series for consistency of variations as function of wavelength and detrending choice. Wrote 
% because dissatisfied with splinelf.m
%
%*** INPUT
%
% User prompted to:
%
% Choose stz-style .mat file (produced by rwlinp.m) or .wk1 spreadsheet for input
% Click on an input .mat file holding the ring-width series (file produced by rwlinp.m)
% Choose either ratio or difference detrending
% Specify number of standard errors for error bars
%
%
%*** OUTPUT 
%
%  Fig 1: Ratio of variance of smoothed index to variance of annual index as function of smoothing, for the three 
%   choices of detrending
%
%
%  For any specified choice of detrending option and smoothing of annual index:
% 
%   Fig 1: Smoothed chronology index, with CI and color-patching of significant departures
%   Fig 2: ure summarizing dependence of between-series agreeement of low freq component on length of detrending spline
%   Figure summarizing variance components due to detrending and smoothing
%   Text variable summarizing common signal in individual seris
% 
%
%*** NOTES
%
% Detrending choices:  
%   Horizontal line at the mean ring width
%   MOdified neg exp; straight line negative slope, or HL at mean if slope positive
%       If ring-width vs time correlation is +, use HL
%       Else try modified exp:  if parameters OK, accept; other wise use SL, neg slope
%   Negative slope combo (NSC):  spline with properties:
%       -everywhere non-increasing
%       -second deriv positive in second half
%       -wavelength with 0.5 amp. freq response no shorter than Ncrit years
%       * OR horizontal line at mean, if above spline not feasible
%   Spline of specified wavelength (by default, the length of the longest series
%
% Smoothing splines.  Ten different gaussian filters are used to smooth the annual series.  The filters have wavelength 
% by default evenly spaced between the length of the longest series and 30 years
%
% 

close all;
clear;
clc;

twarn='Warnings on interpretation of output';
%----   GET THE INPUT RING WIDTHS

kmen1=menu('Choose Format of Input Ring Widths','.mat file via rwlinp.m','Other(not yet implemented)');
nmax=0; % initialize length of longest series
disp('Inputting the ring widths');
if kmen1==1;
    [file1,path1]=uigetfile('*.mat','Input file with ring widths');
    pf1=[path1 file1];
    eval(['load ' pf1 ' X nms yrs;']);
    yrgo = min(yrs(:,1));
    yrsp = max(yrs(:,2));
    yr=(yrgo:yrsp)';
    nyr = length(yr);
    nser = size(nms,1);
    Y=repmat(NaN,nyr,nser);
    % Fill time series matrix
    for n=1:nser;
        igo=yrs(n,1)-yrgo+1;
        nn = yrs(n,2)-yrs(n,1)+1; % number of years of data
        nmax=max([nn nmax]);
        isp = igo+nn-1;
        jgo=yrs(n,3);
        jsp = jgo + nn-1;
        Y(igo:isp,n)=X(jgo:jsp);
    end;
    
    % Rename and clean up
    X=Y;
    nms = cellstr(nms);
    yrs1 = yrs(:,1:2);
    clear yrs  Y n igo isp jgo jsp nn;
    
    % Convert X to mm
    X=X/100;
    LNaN = isnan(X);
    
elseif kmen1==2;
    error('Not yet coded');
    
end;


% Compute the median and minimum series length
nmed=round(median(sum(~isnan(X))));
nmin=min(sum(~isnan(X)));

% Prompt for spline wavelengths
prompt={'Enter min. allowable wavelength for NSC detrending (yr):','Enter wavelength for pure-spline detrending',...
        'Enter min a and max gaussian wavelengths for smoothing'};
def={num2str(nmed),num2str(nmed),num2str([nmed 20])};

% Note the 20-yr minimum gaussian filter. Freq resp of this filter drops to 0.1 at wavelength of 10 yr.
% Thus, the 20-yr gaussian filter essentially smoothes out variations shorter than decadal


dlgTitle='Spline Wavelengths (wavelength at which ampl. of freq response is 0.5 ';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
pnsc=str2num(answer{1});
pspl=str2num(answer{2});
pspair=str2num(answer{3});
ps=linspace(pspair(1),pspair(2),5); % smooth by 5 gaussian filters with wavelengths evenly spaced

% Set some pointers
Lpos=logical(zeros(nser,1)); % allocate for flag for positive linear trend of rw vs time

ntest=nser;

%------ GET TREND LINES 

disp('Estimating growth trend');

% Store growth curves in G{}, raw indice in Y{};  G{j}: j=1:4 -->HL,NE,NSC,SPL

%--- Horizontal line G{1}

disp('  Horizontal Line');
xmn = nanmean(X);
D=repmat(xmn,nyr,1);
D(LNaN)=NaN;
if any(any(D<=0));
    error('Some growth curve goes negative in the HL section');
end;
G{1}=D;
clear D xmn;


%--- Negative exponential G{2}

disp('  Negative Exponential');
D=repmat(NaN,nyr,nser);
for n =1:ntest; % loop over series
    x =X(:,n);
    L1 = ~isnan(x);
    x1 =x(L1);
    nx = length(x1);
    txtyrs =[' (' int2str(nx) ' yrs)'];
    t =(1:nx)';
    r = corrcoef(x1,t);
    r=r(1,2);
    if r>=0;
        g=repmat(mean(x1),nx,1);
        if any(g<=0); 
            error('Somehow got growth curve <=0 under HL option in the NE section');
        end;
        Lpos(n)=1; % flag for positive linear trend with time
        %         plot(t,x1,t,g);
        %         title([nms{n} ' by HL']);
        %         pause(2);
    else;
        [g,k,a,b] = cfnegx1(t,t,x1);
        if ~isempty(g);
            if any(g<=0);
                error(['Neg growth curve for series ' int2str(n) ', NE fit']);
            end;
            %                         plot(t,x1,t,g);
            %                         title([nms{n} ' by mod NE']);
            %                         pause(2);
        else; % must fit neg slope straight line
            F=[ones(nx,1) t];
            b=F\x1;
            g=F*b;
            if any(g<=0); % growth curve goes negative
                % Fit increasingly more flexible splines till g not negative
                delta=0.2; % see splpos.m
                xcrit=min(x);
                [g,w,p,eflag]=splpos(x1,xcrit,t,t,delta);
                if eflag==0;
                    txt1=['Series ' nms{n} txtyrs ': SL gave neg growth curve; switched to ' num2str(w) '-yr spline'];
                    twarn=char(twarn,txt1);
                elseif eflag==1;
                    error( [' Series ' nms{n} txtyrs ': Spline so flexible it blow up']);
                end;
            end;
            %             plot(t,x1,t,g);
            %             title([nms{n} ' by SL']);
            %             pause(2);
        end;
    end;
    D(L1,n)=g;
    
end;
D(LNaN)=NaN;
if any(any(D<=0));
    error('Some growth curve goes negative in the NE section');
end;
G{2}=D;
clear D xmn g b F r t x1 L1 nx;



%--- Negative Slope Combo G{3}

disp('  Negative Slope Combo');
D=repmat(NaN,nyr,nser);

for n =1:ntest; % loop over series
    x =X(:,n);
    L1 = ~isnan(x);
    x1 =x(L1);
    
    nx = length(x1);
    ihalf=round(nx/2);
    
    t =(1:nx)';
    
    % First, using flag from NE method, check that not increasing or zero trend with time.  If such a trend,
    % go right to HL detrend. 
    if Lpos(n); % if positive slope of trend line
        g=repmat(mean(x1),nx,1);
        if any(g<=0); 
            error('Somehow got growth curve <=0 under HL option in the NSC section');
        end;
        D(L1,n)=g;
        %         plot(t,x1,t,g);
        %         title(['Series # ' int2str(n) ': ' nms{n} ' by NSC-->HL']);
        %         pause(2);
    else;
        % Fit spline of maximum acceptable flexibility
        amp=0.5;
        ppnsc=splinep(pnsc,amp);
        g = (cfspl(ppnsc,t,t,x1))';
        d =  diff(g);  %first difference of smoothed curve 
        d2=diff(d); % second difference of smoothed curve
        LL1 = (d <=0.0); % will be 1 if curve never  increasing with time
        LL2=d2(ihalf:(nx-2))>-1e-5; % 1 if change in slope during second half of 
        % positive or imperceptibly negative
        if all(LL1) & all(LL2) & ~any(g<=0);  % curve non-increasing and slope becoming less steep with
            % Slope and change of slope OK, and all growth-curve vales >0. Use it
            D(L1,n)=g;
            %             plot(t,x1,t,g);
            %             title(['Series # ' int2str(n) ': ' nms{n} ' by NSC-->' num2str(pnsc) '-yr spline (max allowable flex)']);
            %             pause(2);
        else;
            % Iteratively find acceptable spline
            [ppnsc,per50] = monotspl(t,t,x1,nx,2);
            g = (cfspl(ppnsc,t,t,x1))'; % Compute spline
            if per50>2*nmax; % if spline wavelength longer than twice series length, use straight line
                F=[ones(nx,1) t];
                b=F\x1;
                g=F*b;
            end;
            if any(g<=0); % growth curve goes negative
                % Fit increasingly more flexible splines till g not negative
                delta=0.2; % see splpos.m
                xcrit=min(x);
                [g,w,p,eflag]=splpos(x1,xcrit,t,t,delta);
                if eflag==0;
                    txt1=['Series ' nms{n} txtyrs ': NSC gave neg growth curve; switched to ' num2str(w) '-yr spline'];
                    twarn=char(twarn,txt1);
                elseif eflag==1;
                    error( [' Series ' nms{n} txtyrs ': Spline so flexible it blow up']);
                end;
                D(L1,n)=g;
                %                 plot(t,x1,t,g);
                %                 title(['Series # ' int2str(n) ': ' nms{n} ' by NSC-->SL']);
                %                 pause(2);
                %                 
                
            else;
                D(L1,n)=g;
                %                 plot(t,x1,t,g);
                %                 title(['Series # ' int2str(n) ': ' nms{n} ' by NSC-->' num2str(per50) '-yr spline ']);
                %                 pause(2);
                
            end;
            
            
        end;
        
    end; % if pos slope
end; % loop over series
if any(any(D<=0));
    error('Some growth curve goes negative in the Neg Slope Combo section');
end;
G{3}=D;



%--- Spline of specified wavelenth G{4}

disp(['  Spline of wavelength ' num2str(pspl) ' yr']);
D=repmat(NaN,nyr,nser);
for n =1:ntest; % loop over series
    x =X(:,n);
    L1 = ~isnan(x);
    x1 =x(L1);
    nx = length(x1);
    t =(1:nx)';
    amp=0.5;
    ppspl=splinep(pspl,amp);
    g = (cfspl(ppspl,t,t,x1))';
    
    if any(g<=0); % growth curve goes negative
        % Fit increasingly more flexible splines till g not negative
        delta=0.2; % see splpos.m
        xcrit=min(x);
        [g,w,p,eflag]=splpos(x1,xcrit,t,t,delta);
        if eflag==0;
            txt1=['Series ' nms{n} txtyrs ': Spline of spec waveL gave neg growth curve; switched to ' num2str(w) '-yr spline'];
            twarn=char(twarn,txt1);
        elseif eflag==1;
            error( [' Series ' nms{n} txtyrs ': Spline so flexible it blow up']);
        end;
    end;
        
    D(L1,n)=g;
    
    %     plot(t,x1,t,g);
    %     title(['Series # ' int2str(n) ': ' nms{n} ' by NSC-->' num2str(pspl) '-yr spline ']);
    %     pause(2);
end;
if any(any(D<=0));
    error('Some growth curve goes negative in the Spline of Specified Wavelength section');
end;
G{4}=D;

%*************  COMPUTE CORE INDICES FROM RING WIDTHS AND GROWTH CURVES

disp('Computing the indices from the ring widths and growth curves');
for n = 1:4;  % loop over the detrendings
    Gthis=G{n};
    Ithis =X ./ Gthis;
    I{n}=Ithis;
end;

%*************  GROUP CORE INDICES BY TREE

disp('Associating cores with trees with treefind.m');
[treenms,itree]=treefind(char(nms));
ntree=length(treenms);
% treenms is col-cell of string names of trees; itree is col vector associating cores to those trees


%************   AVERAGE CORE INDICES INTO TREE INDICES

disp('Averaging core indices into tree indices');
for n =1:4; % loop over detrendings
    Ithis = I{n};  % matrix of core indices
    Tthis = repmat(NaN,nyr,ntree);
    for m = 1:ntree; % loop over trees
        L1 = find(itree==m);
        ncore=sum(L1);
        D=Ithis(:,L1); 
        if ncore==1;
            Tthis(:,m)=D;
        else;
            Tthis(:,m)=(nanmean(D'))';
        end;
    end;
    T{n}=Tthis;
end;


disp('here');    



%************   COMPUTE GAUSSIAN-SMOOTHED TREE INDICES (5 VERSIONS) AND SAMPLES SIZES (# TREES EACH YEAR)


disp('Computing guassian smoothed indices for combinations of detrending and smoothing');


S=cell(4,5); % to hold the smoothed indices
for n =1:4; % loop over detrending
    Tthis=T{n}; % tree indices 
    for m =1:5 % loop over smoothing filters
        pthis=ps(m);
        [Y,U,J]=gausmtx(Tthis,pthis); % Y is the mtx of smoothed indices
        S{n,m}=Y; % store smoothed indices
    end;
end;

% Trial plot the smoothed series for series #3, greatest smoothed, HL vs pure spline detrending
% figure(1);
% uHL = S{1,1}(:,3);
% uSPL = S{4,1}(:,3);
% plot(yr,uHL,yr,uSPL);
% legend('HL detrending',[num2str(pspl) '-yr detrending']);
% title(['Series smoothed by ' num2str(ps(1)) '-yr spline']);
% 

clear Y U J Tthis;


%************   COMPUTE CHRONOLOGY INDEX, SAMPLE SIZE,  AND CI FOR EACH SMOOTHING
disp('Computing average (over trees) smoothed indices, and standard error of mean');
Z = cell(4,5); % to store site-average smoothed chrons
N = cell(4,5); % to store sample size (# trees)
E = cell(4,5); % to store std error of the mean of smoothed indices
SD=cell(4,5); % to store standard deviation 

for m=1:5; % loop over gaussian filters
    for n=1:4; % loop over detrending
        Sthis = S{n,m}; % matrix of smoothed tree indices
        L = ~isnan(Sthis);
        if size(Sthis,2)>1;
            N{n,m}= (sum(L'))'; % col vector of number of trees
            Z{n,m}=    (nanmean(Sthis'))';
            SD{n,m}=    (nanstd(Sthis'))';
        else;
            N{n,m}=numeric(L);
            Z{n,m}= Sthis;
            SD{n,m}=NaN;
        end;
        H=N{n,m};
        Lblow = H==0 | H==1;
        H(Lblow)=NaN;
        E{n,m}=  SD{n,m} ./  sqrt(H);
        clear Lblow;
        
        
    end;
end;



%************  COMPUTE STATISTICS SUMMARIZING EMERGENCE OF 'SIGNAL' FROM NOISE AS SMOOTHING LESSENED
%-max & min index value
%-number of peaks and troughs
%-number of those signif different from 1.0 at  ?SE

%*************  GRAPHICAL OUTPUT


% Build a matrix of choices
dt={'Horizontal Line','Neg Exp','Neg Slope Combo',[num2str(pspl) '-Year spline']};
ds=cellstr(num2str(round(ps')));
dc=cell(4,5);
for n=1:4; % over detrending
    ddt = dt{n};
    for m=1:5; % over smoothing
        dds = ds{m};
        dc{n,m}=['DETREND=' ddt ';  SMOOTH= ' dds  '-yr gaussian'];
    end;
end;
dc=dc(:);        
dc{21}='Back to Previous Menu';

nseq=(1:20)';
jseq=1+floor((nseq-1)/4);
iseq=nseq-(jseq-1)*4;


kwh6=1;
while kwh6;
    kmen1 = menu('Choose','Enter Menus','Quit');
    if kmen1==1;
        kmen2=menu('Choose one',...
            'Time series plot of smoothed chronology and error bars',...
            'Time series plots of families of growth curves');
        if kmen2==1; % 'Time series plot of smoothed chronology and error bars'
            kwh7=1;
            while kwh7;
                kmen3=menu('Choose',dc);
                if kmen3==21;
                    kwh7=0;
                else;
                    i=iseq(kmen3);
                    j=jseq(kmen3);
                    z1=Z{i,j}; % matriz of smoothed indices, detrending i smoothing j
                    e1=E{i,j};
                    figure(1);
                    L=~isnan(z1);
                    tz2 = yr(L);
                    z2=z1(L);
                    e2=e1(L);
                    hp1=plot(tz2,z2,tz2,z2-2*e2,tz2,z2+2*e2);
                    set(gca,'Xgrid','on');
                    xlims=get(gca,'XLim');
                    hline1=line(xlims,[1 1],'Color',[0 0 0]);
                    set(hp1(2),'Color',[1 0 0]);
                    set(hp1(3),'Color',[1 0 0]);
                    title(dc{kmen3});
                end;
            
            end; % while kwh7;
            
        elseif kmen2==2; % 'Time series plots of families of growth curves'
        end;
        
        
    else;
        kwh6=0;
    end;
end; % while kwh6;



uiwait(msgbox('Check the returned variable twarn for any problems with curve fits','Message','modal'));

disp('here');

%************   COMPUTE SAMPLE SIZE (TREE & CORES) TIME SERIES


%*************  COMPUTE CHRONOLOGY INDICES


%*************  COMPUTE GAUSSIAN SMOOTHED INDICES (5 VERSIONS)





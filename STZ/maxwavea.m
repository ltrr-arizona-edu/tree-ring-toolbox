function maxwave(kquick)
% maxwave: assess maximum wavelength of intercorrelation of tree-ring indices between trees or cores
%
%
%*** INOYT
%
% kquick (1 x 1)i  option of quick mode; ==1 yes, ==0 no
% SUGGEST CHANGES
%
% Option to view median stat of each core with all others -- a color coded and labeled figure
% Finish option for series info -- maybe use above for combo into
% This a big change:  option for trial and error masking of series at start of program.  Ideally,
%  by point and click to series box ; combine with appendable stored mask vector that can eliminate
%  need to re-mask each time start finction.  This mask can be stored in the main .mat matrix. If
%  want, can ignore the mask, or can change it.

% 
%*** NOTES
%
% Objective,  a color figure for each of several target wavelengths (e.g., 200 yr, 300 yr)
% Figure will have:
%   Series number along each axis
%   Correlation between series pair as color coded cell value
%   Color bar key to correlation
%   Title with target wavelength, names of series
% Strategy.
%   Use spline1 for detrending rw by ratio method -->index
%   Use spline2 for smoothing the index to emphasize target wavelength and longer
% Steps
%   Prompt for specs
%   Get ring widths
%   Cull series long enough to show the target wavelength (at least 1 cycle)
%   Compute spline parameter of splines with freq response 0.1 and 0.9 at the target wavelength
%   Detrend series with spline1 
%   Convert series to indices
%   Smooth indices with spline2
%   Cross-correlate series
%   Display map of cross-correlation in one window
%   Display cross-ref list to nms in a text window
%   Display mix of all smoothed indices in one figure
%   Optionally plot smoothed series times series plots in a third window
%
%-- Figure Windows
% 1 text info on program settings
% 2 change of similarity stat with target wavelength
% 3 summary over target wavelengths
% 4 time series plots, smoothed pair -- a selected pair, for current target wavelentgr
% 5 reserved for tsp of indiv series vs robust lfi based on others -- current target wavelength
% 6:(npd+6)  color matrices for target wavelengths
%
% 

% Hard Code
close all;
shading flat;

txtsort={'Sequence number','Length','First year','Alphabetical series id'};
txtstat ={'Correlation','Cross product of departures from 1.0',...
        'Relative difference','Coefficient of departure'};

ngofig = 6; %First figure window for a color matrix


%   PROMPT FOR SPECS

% Load input file
[file1,path1]=uigetfile('*.mat','Infile with ring widths');
pf1=[path1 file1];
eval(['load ' pf1 ' nms X yrs maskwv;']); % load ring widths
if  ~all([exist('nms','var')==1 exist('X','var')==1  exist('yrs','var')==1]);
    error([pf1 ' does not contain X, nms ans yrs']);
end;
nser = size(nms,1); % number of series, before any masking
d = yrs(:,2)-yrs(:,1)+1; % series lengths
nmsall=nms;
yrsall=yrs;
dall=d;


% Handle masking.  You may have previously set elements of maskwv to 0, indicating those series not
% to be used in the analysis
Lmask = logical(ones(nser,1));
% Check for mask of series from analysis
if exist('maskwv','var')==1;
    Lmasknow = maskwv; % store logical mask as Lmasknow
    kignore = questdlg('Ignore series mask');
    switch kignore;
    case 'Yes';
        strtemp='Number of series (none masked) = ';
        % Lmask still logical vector all ones
        % no action
    otherwise; % omit masked series from analysis
        Lmask=Lmasknow;
        yrs = yrs(Lmask,:);
        nms = nms(Lmask,:);
        nmask=sum(Lmask==0);
        strtemp=['Number of series (minus ' int2str(nmask) ' masked) = '];
        nser = size(nms,1);
    end;
else;
    % no action
    strtemp='Number of series (no mask) = '; 
end;
I1= find(Lmask); % relative index of each of the currently ordered and selected series into  nmsall, yrsall, dall 



% Lengths of series, and default length of a "long" series.  This default used in prompt allowing user
% to restrict analysis to long series
d=dall(Lmask); 
dlong = prctile(d,75); % 75 th percentile of series lengths 
dlong = 10*ceil(dlong/10);


strtxt1={['INPUT RING-WIDTH FILE = ' file1],...
        [strtemp num2str(nser)],...
        ['   Shortest = ' num2str(min(d)) ' yr'],...
        ['   Longest = ' num2str(max(d)) ' yr'],...
        ['   Earliest = A.D. ' num2str(min(yrs(:,1)))],...
        ['   Latest = A.D. ' num2str(max(yrs(:,2))) ]};
clear strtemp;



% ----  TARGET WAVELENGTHS

% Set default longest target wavelength == nearest number divisible by 100 shorter than
% second longest series (because will eventually need at least two series as long as each target 
% wavelength

% Longest 
dsort=sort(d);
ndsort=length(dsort);
pd = dsort(ndsort-1); % second longest series
pd=100* floor(pd/100);  
pdstr=num2str(pd,'%4.0f');

% Decrement
if pd>1000
    delta1 = 200;
else;
    delta1=100;
end;
delta1str = num2str(delta1,'%4.0f');

% Shortest wavelength
waveshort = 20; 
waveshortstr = num2str(waveshort,'%4.0f');

% Prompt
prompt={'Longest (yr):','Decrement','Shortest'};
def={pdstr,delta1str,waveshortstr};
dlgTitle='Target Wavelengths (yr)';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
pdstr=answer{1};
pd=str2num(pdstr);
delta1str = answer{2};
delta1 = str2num(delta1str);
waveshortstr = answer{3};
waveshort = str2num(waveshortstr);

if waveshort<20;
    error('Shortest allowable target wavelength 20 yr');
end;
if pd>=100; % if longest target wavelength exceeds 100 yr
    if (length(pd:-delta1:100) + 2)>20;
        error('More than 20 target wavelengths');
    end;
    if waveshort<50;
        pdvec = [pd:-delta1:100  50 waveshort];
        
    elseif waveshort<100;
        pdvec = [pd:-delta1:100  waveshort];
    else; 
        pdvec=[pd:-delta1:waveshort];
    end;
else; % longest target wavelength less than 100 yr
    pdvec=[pd:-delta1:waveshort];
    if length(pdvec)>20;
        error('More than 20 target wavelengths');
    end;
end;
npd=length (pdvec); % number of target wavelengths

strtemp = {'  ',...
        ['TARGET WAVELENGTHS (yr)'],...
        ['   Longest = ' pdstr],...
        ['   Decrement = ' delta1str],...
        ['   Shortest = ' waveshortstr],...
        ' ',...
        ['RESTRICTIONS']};
clear dsort ndsort;
strtxt1 = [strtxt1 strtemp];
clear strtemp dsort ndsort; 



% OPTIONAL RESTRICTION OF ANALYSIS TO LONG SERIES 

klong = questdlg('Restrict analysis to long series');
switch klong;
case 'Yes';
    prompt={'Enter minimum allowable series length (yr):'};
    def={num2str(dlong)};
    dlgTitle='Restrict analysis to series at least this long';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    nlong=answer{1};
    strtemp = {['Use only series with length at least ' num2str(nlong) ' yr']};
    strtxt1 = [strtxt1 strtemp];
    nlong = str2num(nlong);
otherwise;
end; % switch klong



% MINIMUM OVERLAP FOR COMPUTATION SIMILARITY STATISTIC
prompt={'Enter minimum acceptable overlap (yr):'};
def={'50'};
dlgTitle='Minimum Acceptable Overlap for Correlation, etc. (yr)';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
minss=answer{1};
nmin=str2num(minss);
if nmin >pd | nmin<10;
    error(['Minimum overlap must not be less than 10 yr or greater than ' pdstr ' yr']);
end;
strtemp = {['Minimum acceptable overlap at least '  num2str(nmin) ' yr']};
strtxt1 = [strtxt1 strtemp];
    

% Option for correlation matrix vs relative departire statistic
%
% I have decided to compute the RD statistic this way:
%
% MSD1 = mean square departure of smoothed indices for series a and series b
% MSD2 = mean square departure of smoothed indices from normal (==1.0)
%  Could compute a MSD2 value for each series and take average.  I prefer compute the 2N squared departures
%  and compute the mean of these.  Not sure if this is same as harmonic mean of the two individual means or what
% RD = 1 - (MSD1/MSD2)
%
% Reasoning. Want measure of similarity of the two smoothed indices.  Both could be way below normal, have 
% slightly opposing trends for the ovelap, and a correl coef of near -1.0.  The RD avoids this.  If both smoothed series
% indicate very low growth, both way below 1.0, then the RD above is likely to be positive.  
% Also, smoothing invalidates conventional test of signif of correl coeff, so nothing lost there in using RD as measure

khow=menu('Choose Method for Similarity','Correlation','Median cross-product','Relative Difference','Coefficient of Departure');
if khow==1;
    methsim='R';
    strmeth='Correlation';
elseif khow==2;;
    methsim='CP';
    strmeth = 'Cross-Product';   
elseif khow==3;;
    methsim='RD';
    strmeth = 'Relative Difference';
elseif khow==4;
    methsim='CD';
    strmeth='Coefficient of Departure';
end;

strtemp = {' ',  ['SIMILARITY STATISTIC:  ' txtstat{khow}]};
strtxt1 = [strtxt1 strtemp];


% SORTING OF SERIES 

ksort = menu('Choose Sorting Method',...
    'Sequence as in nms in input .mat file',...
    'Length (1==longest)',...
    'Age (1==earliest first year',...
    'Alphabetically by series id');
if ksort==1; % Sequence as in nms in input .mat file
    % no action needed
    strtemp={'Series not re-ordered'};
    Is = (1:nser)';
elseif ksort==2; % Length (1==longest)
    [s,Is]=sort(d);
    s=flipud(s);
    Is=flipud(Is);
    strtemp={'Series sorted by length (1=longest)'};
elseif ksort==3; % Age (1==earliest first year
    [s,Is]=sort(yrs(:,1));
    strtemp={'Series sorted by age (1=earliest)'};
elseif ksort ==4; % Alphabetically by series id
    [s,Is]=sort(cellstr(nms));
    s=char(s);
    strtemp={'Series sorted alphabetically by id'};
end;
if ksort~=1; % if re-sorting data
    nms=nms(Is,:);
    yrs=yrs(Is,:);
    d=d(Is,:);
end;

strtxt1 = [strtxt1 strtemp];




Lmasknew = logical(ones(nser,1));
kwh6=1;
while kwh6;  % while loop for interactive masking
    
    kquest7=questdlg('Run, or re-run with revised mask');
    switch kquest7;
    case 'Yes';
        
        
        % Allocate cells to hold results for each target wavelength
        Scell = cell(1,npd); % allocate to store time series of smoothed time series
        Rcell = cell(1,npd); % correlations between series
        Ncell = cell(1,npd); % Sample size (# yr) correls in Rcell based on
        nmscell = cell(1,npd); %  names of series in Rcell
        tcell=cell(1,npd); % time vector for matrices in Scell
        Icell=cell(1,npd); % row cross-ref of series to original order in nms in input .mat file
        
        
        % LOOP OVER TARGET WAVELENGTHS
        
        for k = 1: npd;
            pd=pdvec(k);
            
            
            %   CULL SERIES LONG ENOUGH TO SHOW TARGET WAVELENGTH (1 CYCLE)
            
            Lcull = d >= pd;
            if strcmp(klong,'Yes');
                Llong = d>=nlong; 
                Lcull = Llong & Lcull & Lmasknew;
            else; 
                Lcull=Lcull & Lmasknew;
            end;
            
            
            nkeep=sum(Lcull);
            if nkeep<2;
                error(['Fewer than 2 series of length ' pdstr ' yr']);
            end;
            Ic = find(Lcull); % relative index
            % nmsall(I1(Is(Ic))) gives names
            
            % Earliest and latest year of culled series
            yrgo = min(yrs(Lcull,1));
            yrsp = max(yrs(Lcull,2));
            yrY = (yrgo:yrsp)'; % year vector for matrix Y
            nyr = length(yrY);
            
            % Culled series into matrix
            Y=repmat(NaN,nyr,nkeep);
            i1 = yrs(Lcull,1)-yrgo+1; % target first row of each series in Y
            i2 = i1+ (yrs(Lcull,2)-yrs(Lcull,1)-1); % last row
            % starting and ending rows of source data in X
            j1 = yrs(Lcull,3);
            j2= j1 + (yrs(Lcull,2)-yrs(Lcull,1)-1);
            
            
            for n = 1:nkeep;
                x=X(j1(n):j2(n));
                Y(i1(n):i2(n),n)=x;
            end;
            
            
            
            %   COMPUTE  spline parameter of splines with freq response 0.1 and 0.9 at the target wavelength
            p1=splinep(pd,.1); % detrending spline, little tracking of variations at target period
            p2=splinep(pd,.90); % smoothing spline: strong tracking
            
            pd50a = spline50(p1); % wavelength at which freq response of spline1 is 0.5
            pd50b = spline50(p2); % wavelength at which freq response of spline1 is 0.5
            if pd50a<0;
                pd50a=-pd50a;
            end;
            if pd50b<0;
                pd50b=-pd50b;
            end;
            str50a = num2str(pd50a);
            str50b = num2str(pd50b);
            
            
            
            %   DETREND  CULLED SERIES WITH SPLINE1
            
            D=Y; % to hold trend lines
            % Compute detrending spline
            for n = 1:nkeep; % loop over series
                yall =Y(:,n);
                Lgood=(~isnan(yall));
                y=yall(Lgood);
                yry=yrY(Lgood);
                g=csaps(yry,y,p1,yry);
                D(i1(n):i2(n),n)=g;
            end;
            
            
            
            %   CONVERT TO INDICES
            
            I = Y ./ D;
            
            
            
            %   SMOOTH INDICES
            
            S=Y; % to hold trend smoothed indices
            
            for n = 1:nkeep; % loop over series
                yall =I(:,n);
                Lgood=(~isnan(yall));
                y=yall(Lgood);
                yry=yrY(Lgood);
                s=csaps(yry,y,p2,yry);
                
                % Stored smoothed index
                S(i1(n):i2(n),n)=s; 
            end;
            
            
            
            % CALL FUNCTION TO COMPUTE INTER-SERIES CORRELATIONS OR RELATIVE DEPARTURE
            
            if strcmp(methsim,'R');
                [R,N]=simptch(S,nmin,1);
            elseif strcmp(methsim,'CP');
                [R,N]=simptch(S,nmin,2);
            elseif strcmp(methsim,'RD');
                [R,N]=simptch(S,nmin,3);
            elseif strcmp(methsim,'CD');
                [R,N]=simptch(S,nmin,4);    
            end;
            
            
            % Make  matrix symmetric
            R1 = rot90(R);
            R1=flipud(R1);
            A=ones(nkeep);
            LA=logical(tril(A,-1));
            R(LA)=R1(LA);
            
            Scell{k}=S;
            Rcell{k}=R;
            Ncell{k}=N; 
            nmscell{k}=nms(Lcull,:);
            tcell{k}=(yrgo:yrsp)';
            Icell{k}=(Ic);
            
            
            % COLOR MATRIX PLOT
            
            if kquick~=1;
                
                figure(ngofig+k-1);
                pdthis=pdvec(k);
                strtit = [strmeth ' for Target Wavelength ' num2str(pdthis) ' Yr'];
                
                t = tcell{k};
                nmthis=nmscell{k};
                nthis=Ncell{k};
                Sthis=Scell{k};
                
                [mR,nR]=size(R);
                % To use pcolor, need first and last row of NaN
                zr = repmat(NaN,1,nR);
                R=[R; zr];
                zc= repmat(NaN,mR+1,1);
                R=[R zc];
                
                pcolor(R);
                
                colord=get(gca,'ColorOrder');
                set(gca,'XTick',(1:mR)+0.5);
                set(gca,'XTickLabel',num2str((1:mR)'));
                set(gca,'YTick',(1:nR)+0.5);
                set(gca,'YTickLabel',num2str((1:nR)'));
                
                b1temp = abs(nanmin(nanmin(R)));
                b2temp = abs(nanmax(nanmax(R)));
                b3temp=max([b1temp b2temp]);
                climit = 1.2*[-b3temp b3temp];
                clear b1temp b2temp b3temp;
                set(gca,'CLim',climit);
                colorbar;
                xlabel('x Series');
                ylabel('y Series');
                
                title(strtit);
                
            end;
            
            disp(['Finished storing results for target wavelength ' num2str(pd) ' yr']);
            
            
            
            
        end; % for k = 1: npd;
        
        
        
        
        % COMPUTE AND STORE CORREL CHANGE WITH WAVELENGTH
        
        % recall that npd is number of periods
        rset=cell(npd,1); % to hold nonredundant corrs
        
        minr = 0; % Initialize y limits for summary similrity plots
        maxr=0;
        for n = 1:npd; % loop over target wavelengths
            R=Rcell{n}; % correl matrix
            [mR,nR]=size(R);
            
            %     U = ones(mR);
            %     U=logical(triu(U,1));
            %     Rvec = R(U);
            %     Rvec=Rvec(~isnan(Rvec));
            Rvec=nanmedian(R);
            Rvec=Rvec(~isnan(Rvec));
            minr=min([Rvec minr]);
            maxr=max([Rvec maxr]);
            rset{n}=Rvec;
        end;
        if strcmp(methsim,'R'); % If correlation coef, y limit is -1 to 1
            minr=-1.0;
            maxr=1.0;
        end;
        
        clear R mR nR U Rvec ;
        
        
        
        % CHANGE OF COREL WITH WAVELENGTH PLOT  
        
        figure(2);
        hax1 = axes;
        set(gca,'XLim',[0 (npd+1)]);
        set(gca,'YLim',[minr maxr]);
        hold on;
        
        
        
        for n = 1:npd; % over wavelengths
            rvec = rset{n}; % correlations
            nvec=length(rvec); % sample size
            plot(repmat(n,nvec,1),rvec,'.',n,median(rvec),'o');
            
        end;
        if ~exist('colord')==1;
            colord=get(gca,'ColorOrder');
        end;
        hline3=line([0 (npd+1)],[0 0]);
        set(hline3,'Color',colord(3,:));
        hold off;
        xlabel('Target Wavelength (yr)');
        ylabel(strmeth);
        title(['Median ' strmeth ' at Target Wavelengths']);
        xtL = cellstr(num2str([NaN pdvec NaN]'));
        xtL{1}='';
        xtL{npd+2}='';
        set(gca,'XTick',[0:1:(npd+2)]);
        set(gca,'XTickLabel',xtL);
        
        
        %--  SUMMARY FOR SERIES OVER TARGET WAVELENTHS 
        
        % Allocate 
        W = repmat(NaN,npd,nser); % to hold statistic for color plot
        
        for n = 1:npd; % loop over target wavelengths
            q = rset{n};
            icell = Icell{n};
            nmcell = nmscell{n};
            W(n,icell)=q;
            % Depending on screening criteria, some cols of W may be all-NaN.  Want to avoid that,
            % but if delete those cols, must keep track of indices of remaining cols into nms
            
        end;
        jgood = ~all(isnan(W)); 
        I4 = find(jgood);  % note that nms(I4,:) will recover the series names
        ngood = length(I4);
        W=W(:,I4);
        [mW,nW]=size(W);
        Wnms = nms(I4,:)
        figure(3);
        if kquick~=1;
            zw = repmat(NaN,1,nW);
            W=[W; zw];
            zc= repmat(NaN,mW+1,1);
            W=[W zc];
            
            pcolor(W);
            set(gca,'XTick',(1:nW)+0.5);
            set(gca,'XTickLabel',num2str((1:nW)'));
            set(gca,'YTick',(1:mW)+0.5);
            set(gca,'YTickLabel',num2str((1:mW)'));
            
            b1temp = abs(nanmin(nanmin(W)));
            b2temp = abs(nanmax(nanmax(W)));
            b3temp=max([b1temp b2temp]);
            climit = 1.2*[-b3temp b3temp];
            clear b1temp b2temp b3temp;
            set(gca,'CLim',climit);
            colorbar;
            xlabel('');
            ylabel('Target Wavelength (yr)');
            set(gca,'Position',[.13 .31 .6626 .51]);
            set(gca,'YTickLabel',num2str(pdvec'));
            set(gca,'XTickLabel',[]);
            % Add series names along x axis
            xpt = get(gca,'XTick');
            ypt = get(gca,'YLim');
            ypt = ypt(1);
            for nn = 1:ngood;
                wname = Wnms(nn,:);
                xpt1 = xpt(nn);
                text(xpt1,ypt,wname,'Rotation',60,'HorizontalAlignment','right');
            end;
            
        end;
        % --- TEXT WINDOW
        
        figure(1);
        text(0.1,0.95,strtxt1,'VerticalAlignment','Top');
        set(gca,'Visible','off');
        title('PROGRAM SETTINGS');
        
        
                
        kwh1=1;
        
        while kwh1==1; 
            kmen5='No';  % initialize for no additional masking yet done
            figure(1);
            kmen1=menu('Choose One',...
                'Color mapped similarity of each series with all others for a target bandpass wavelength',...
                'Change in similarity statistic as function of target bandpass wavelength',...
                'Color mapped summary of similarity for each series',...
                'Quit');
            if kmen1==1; % color maps, each series vs all others
                
                kwh1=1;
                kwh2=1;
                while kwh2==1;
                    
                    men2chc = cellstr(num2str(pdvec'));
                    men2chc{npd+1}='Quit';
                    kmen2 = menu('Choose Target Wavelength (yr):',men2chc);
                    
                    if kmen2==npd+1;
                        kwh2=0;
                    else;
                        
                        figure(ngofig+kmen2-1); % bring color matrix as current figure
                        pdthis=pdvec(kmen2);
                        strtit = [strmeth ' for Target Wavelength ' num2str(pdthis) ' Yr'];
                        
                        
                        R=Rcell{kmen2};
                        t = tcell{kmen2};
                        nmthis=nmscell{kmen2};
                        nthis=Ncell{kmen2};
                        Sthis=Scell{kmen2};
                        icull=Icell{kmen2};
                        
                        [mR,nR]=size(R);
                        
                        
                        kwh3=1;
                        while kwh3==1;
                            
                            kquest1=questdlg('Get series-pair details by clicking on cell');
                            switch kquest1;
                            case 'Yes';
                                figure(ngofig+kmen2-1);
                                [x1,y1]=ginput(1);
                                % Calculate which series those are
                                x1=round(x1-0.5);
                                y1=round(y1-0.5);
                                
                                namex = nmthis(x1,:);
                                namey= nmthis(y1,:);
                                sx = Sthis(:,x1); % smoothed series
                                irefnear = (icull(x1)); % reference to nms row in sorted, stored-mask series
                                ireffar = Is(irefnear); % reference to nms row in nmsall
                                sy = Sthis(:,y1);
                                Lx = ~isnan(sx);
                                Ly = ~isnan(sy);
                                yron = max([min(t(Lx))  min(t(Ly))]);
                                yroff = min([max(t(Lx))  max(t(Ly))]);
                                if yron>yroff;
                                    uiwait(msgbox([namex ' does not overlap ' namey],'Title','modal'));
                                    kwh4=0;
                                elseif (yroff-yron)<nmin;
                                    uiwait(msgbox([namex ' has too small an overlap with ' namey],'Title','modal'));
                                    kwh4=0;
                                else;
                                    kwh4=1;
                                end;
                                
                                tyr = (yron:yroff)';
                                Lt = t>=yron & t<=yroff;
                                sx = sx(Lt);
                                sy = sy(Lt);
                                
                                
                                while kwh4==1;
                                    kmen4=menu('Choose','Time series plots of smoothed indices',...
                                        'Series information',...
                                        'Quit');
                                    if kmen4==1;
                                        figure(4);
                                        
                                        % Compute the selected similarity statistic for display in title with correl coef
                                        % First, correl coef
                                        rtemp=corrcoef(sx,sy);
                                        rtemp=rtemp(1,2);
                                        rtempstr= ['r = ' num2str(rtemp,'%5.2f')];
                                        % Next, the selected statistic
                                        switch methsim;
                                        case 'R';
                                            ss=[];
                                        case 'CP';
                                            ss = 100*median ((sx-1.0) .* (sy-1.0)) 
                                        case 'RD';
                                            msa = mean((sx - sy).^2);
                                            msb = mean([(sx-1.0) ; (sy-1.0)] .^2);
                                            ss = 1 - (msa/msb);
                                        case 'CD';
                                            ss= mean((sx-1.0) .* (sy-1.0))/(std(sy)*std(sx));
                                        end;
                                        ssstr = [methsim ' = ' num2str(ss)];
                                        
                                        
                                        hp2=plot(tyr,(sx),tyr,(sy));
                                        set(hp2(2),'LineStyle','--');
                                        
                                        hline2=line([min(tyr) max(tyr)],[1 1 ]);
                                        set(hline2,'Color',colord(3,:));
                                        set(hp2(1),'LineWidth',2);
                                        set(hp2(2),'LineWidth',2);
                                        xlabel('Year');
                                        ylabel('Smoothed Index');
                                        legend(namex,namey);
                                        grid;
                                        switch methsim;
                                        case 'R';
                                            title(['Target Wavelength =  ' num2str(pdthis) ' yr;   '  rtempstr  ]);
                                        otherwise;
                                            title(['Target Wavelength =  ' num2str(pdthis) ' yr;   '  rtempstr '; ' ssstr ]);
                                        end;
                                        
                                        
                                        
                                        
                                        
                                    elseif kmen4==2;
                                    elseif kmen4==3;
                                        kwh4=0;
                                        figure(2);
                                    end;
                                end; % kwh4==1;
                                
                                
                            otherwise;
                                kwh3=0;
                                
                            end;
                        end; % kwh3==1;
                        
                    end;
                    
                end; % kwh2==1;
                
                
            elseif kmen1==2;  % change in correlation
                figure(2);
                uiwait(msgbox('Click OK to continue','Message','modal'));
                
                kwh1=1;
            elseif kmen1==3; % summary similarity for series/wavelength
                figure(3);
                strtit3={'MEDIAN SIMILARITY FOR EACH SERIES',...
                        ['Statistic = ' txtstat{khow}]};
                
                title(strtit3);
                kwh8=1;
                while kwh8;
                    kmen8=menu('Choose',...
                        'Back up one menu',...
                        'Show series info (click using crosshair)',...
                        'Mask series from subsequent analysis');
                    if kmen8==1; % up one menu
                        kwh8=0;
                    elseif kmen8==2; % Show series info
                        [x1,y1]=ginput(1);
                        x1=round(x1-0.5);
                        y1=round(y1-0.5);
                        irefnear = (I4(x1)); % reference to nms row in sorted, stored-mask series
                        ireffar = Is(irefnear); % reference to nms row in nmsall
                        
                        wname = Wnms(x1,:);
                        wperiod = ['Period = ' num2str(yrs(irefnear,1:2))];
                        uiwait(msgbox({wname,wperiod},'Message','modal'));
                        
                    elseif kmen8==3; 
                        kmen5='No';
                        [x1,y1]=ginput(1);
                        x1=round(x1-0.5);
                        y1=round(y1-0.5);
                        irefnear = (I4(x1)); % reference to nms row in sorted, stored-mask series
                        ireffar = Is(irefnear); % reference to nms row in nmsall
                        wname = Wnms(x1,:);
                        kwh9=1;
                        while kwh9==1;
                         
                            kmen9=menu('Choose',...
                                'Mask it',...
                                'Nah,do not mask it');
                            if kmen9==1; % desire to mask
                                Lmasknew(irefnear)=0;
                                maskwv(ireffar)=0;
                                kmen5='Yes';
                            else;
                            end;
                        end;
                      
                    end;
                end;
                kwh1=1;
                    
                        
            elseif kmen1==4;
                kwh1=0;
            end;
            %         if strcmp(kmen5,'Yes');
            %             kwh1=1;
            %         else;
            %             kwh1=0
            %         end;
            
            
        end; % kwh1==1; 
        if strcmp(kmen5,'Yes');
            kwh6=1;
        else;
            kwh6=0
        end;
        
    case 'No';
        kwh6=0;
    otherwise;
    end;
    
    
end;  % while kwh6==1  % masking loop


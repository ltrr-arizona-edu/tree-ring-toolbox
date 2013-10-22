function [Y,f]=sprdloc(X,lenseg,mlap,txtin,fprev,kopt)
% sprdloc: matched power transformation guided by slope of spread-vs-level plot
% [Y,f]=sprdloc(X,lenseg,mlap,txtin,fprev,kopt);
% Last revised 1-29-01
%
% A power transformation can be used to adjust or stabilize the variance as a function of level, or location, in a time series. 
% Tree-ring standardization usually incorporates such adjustment, either by computing indices as ratios of ring width to a fitted
% line or by log-transforming ringwidth before computing indices as the difference of ringwidth and a fitted line.
% Sprdloc.m can be used to check the suitability of various power transformations for stabilizing variance.  Sprdloc.m allows 
% automated or interactive graphical manual selection of the appropriate power using as a guide the slope of the 
% spread-vs-level plot.  Automated modes allow user to specify log transform, have the power computed directly from
% the slope of the spread-vs-level plot, or have the power selected from a suite of common power transforms using the
% slope of the spread-vs-level plot as a guide.  This last mode promotes conservative transformation by preferring the most 
% conservative power transformation compatible with the slope of the spread-vs-level plot and its uncertainty.
% As a final step, the transformed series is matched, or linearly re-scaled, so that the units resemble those of the original
% series.  Matching is performed such that the minimum and median of the transformed series equal the minimum and median of the 
% original series.  This scaling is convenient for ring-width interpretation because locally absent rings are zero in both the 
% original and transformed data, and the scale of transformed ring width is similar to that of measured ring width.
%
%*** INPUT ARGUMENTS *************************
%
% X (mX x 2)r  time series matrix, with year as column 1 and data as column 2. Data must be non-negative (see notes)
% lenseg (1 x 1)i length of averaging window used to define batches (e.g., 30 yr) (see notes)
% mlap (1 x 1)i  overlap of segments used to define batches (e.g., 10 yr)
% txtin{1 x 3} text input
%   {1}(1 x ?)s id of time series (e.g., FTP01A)
%   {2}(1 X ?)s variable (e.g., Ring width)
%   {3}(1 x ?)s units (e.g., mm x 100)
% fprev{}(1x6) information on previous matched transformation ([] is this is first attempt) (see notes)
%   {1} c (1 x 1)r  shift parameter applied to input series x before transformation (x+c is transformed)
%   {2} p (1 x 1)r  power of transformation applied to the shifted input (T(x+c)=(x+c)^p )
%   {3} a (1 x 1r  first parameter, a,  applied in matching:  y=a + bT(x+c), where T(x+c) is the transformation
%   {4} b (1 x 1)r second parameterm b, in matching: y = a + bT(x)
%   {5} eqnstr(1 x ?)s string variable for transformation (e.g, y=x^{.5} for square root)
%   {6} khow (1 x 1)i how the transformation power selected (see output for key)
% kopt (1 x 1)i options
%   kopt(1) fit mode
%       ==1 log-transform (base 10)
%       ==2 commonly used power compatible with slope of spread vs level (see notes)
%       ==3 power computed from slope of spread vs level
%       ==4 interactive fitting
%
%*** OUTPUT ARGUMENTS ************************
%
% Y (mY x 2)r transformed time series, with year as col 1
% f{}(1x5) information on selected matched transformation 
%   {1} c (1 x 1)r  shift parameter applied to input series x before transformation (x+c is transformed)
%   {2} p (1 x 1)r  power of transformation applied to the shifted input (T(x+c)=(x+c)^p )
%   {3} a (1 x 1r  first parameter, a,  applied in matching:  y=a + bT(x+c), where T(x+c) is the transformation
%   {4} b (1 x 1)r second parameterm b, in matching: y = a + bT(x)
%   {5} eqnstr(1 x ?)s string variable for transformation (e.g, y=x^{.5} for square root)
%   {6} khow (1 x 1)i how the transformation power selected (see output for key)
%       ==1 automated: log (base 10)
%       ==2 automatic: commonly used power (one of eight possible) compatible with slope of spread vs level
%       ==3 automatic:  power computed as 1-slope, where slope is slope of of log(spread) vs log(median)
%       ==4 Manual, accepting p given by slope of spread vs level
%       ==5 Manual, specified p overrides automatic p
%
%*** REFERENCES 
%
% Hoaglin, D.C., Mosteller, F., and Tukey, J.W., 1983, Understanding Robust and Exploratory Data Analysis: 
% John Wiley & Sons, Inc.,  New York, 447p.
%
%*** UW FUNCTIONS CALLED
%
% hinge -- divides time series into batches and computes hinges (minimim,lower quartile, median, upper quartile, max)
% shiftup1 -- utility function to shift if needed to avoid zero values in series to be transformed
% subfn1 -- subfunction that computes "matched" transformation
%
%*** TOOLBOXES NEEDED-- NONE
%
%*** NOTES ***********************
%
% The matched transformation is given by
%   y = a + b*T(x+c),
% where x is the original series, c is a shift parameter, T(x) is the log transform or power transform of x, 
% and a and b are coefficients computed such that the median and minimum of y are equal to the median and minimum of x.
%
% Idea is to compute power of tranformation from the slope of the spread-vs-level plot, a plot of log(Q) vs log(M),
% where Q is the interquartile range and M is the median for various "batches" of observations. 
% The batches are overlapping segments of the time series.  The interquartile range is the difference between 
% the upper and lower "fourths" as defined by Hoaglin et al. (1983). The slope of the spread vs level plot is directly 
% related to the power of the appropriate transformation by p = 1-b , where p is the power and b is the slope. If p is
% approximately zero, the log transformation is indicated.  
%
% The slope may have a wide confidence interval, making the choice of power uncertain.  The automated fitting option
% under kopt(1)==2 handles this problem by considering the approximate 95% confidence interval of the slope and preferring
% the most conservative of several transformations compatible with the uncertain slope. For example, if the 95% CI of
% the slope includes zero, the null transformation (no transform) is selected.
%
% Input time series x must be non-negative. x may have zero values (e.g., ringwidth measurements zero), but if
% so, x may be shifted before power transformation.  To shift, a small increment is added to each data value.  The 
% increment equals the smallest nonzero value in the original time series.  The output argument c is this increment.
% Series with zero values are shifted in this way only if the selected transformation is a log-transform or a 
% power transformation with a negative exponent. x is permitted to have NaNs imbedded.  Thus, ring-width series with 
% some internal undatable segment can be handled.
%
% Reversal of sign of transformed series if power negative.  If the selected power is negative (e.g., -1 for an inverse 
% transform) the sense of anomalies is reversed from that of the original series. For example, the largest value in the original 
% series x becomes the smallest transformed value 1/x.  To return the transformed data to the sense of the original, the
% transformed series is reversed in sign if the power is negative.  Thus, for negative power p, the transformation is
% y = a + b*(-x^p),  where ^ denotes exponentiation  
%
% lenseg and mlap.  These input arguments control how "batches" of data are defined.  Batches in sprdloc.m are groups of 
% successive observations in a time series.  Sprdloc.m uses a moving-average window to define the batches. The window has
% length lenseg years.  The batches may overlap, as specified by the parameter mlap (yrs).  For example, settings of 
% lenseg=30 and mlap=10 specify that batches are 30-year groups of data overlapped by 10 years.  Overlapping batches may help
% identify time trends in spread and location, but has the undesirable effect of invalidating the confidence interval on the
% computed splope of the spread-vs-level plot.  The problem is that the obsevations (batches) are not independent -- adjacent 
% batches have some of the same data.  The confidence intervals should therefore be wider than computed. The confidence bands
% are nonetheless useful as a rough guide to avoiding power transformations unjustified by the strength of the spread vs 
% level relationship.
%
% Common power transformations.  Option kopt(1)=2 specifies that sprdloc.m automatically chooses from among several available 
% common power transformations.  Most candidates are members of Tukey's "ladder of powers" listed in Table 3-10 
% of Hoaglin et al. (1983, p. 82).  The candidate transformations are 1) cube, 2) square, 3) none, 4) p=2/3, 
% 5) square root, 6) p=1/3,  7) log10, 8) p=-1/3,  9) p=-1/2 (reciprical root),  10) p=-2/3,  11) reciprical, and 
% 12) reciprical square.   Most conservative is power=1, or no transformation.  Powers
% closer to 1 are considered more conservative than powers further from 1.  The most conservative power within the 99% 
% confidence interval of the computed slope is selected. If none of the common powers are in the confidence interval, the
% closest more conservative power is selected.
%
% From experimenting with randomly generated data, I found that strictly accepting the power computed from the slope of the
% spread-vs-level plot is not a good idea.  Too often transformation is indicated where the data are generated from
% a distribution with a specified mean and variance.  The automated mode (kopt(1)=2) helps reduce this problem by the 
% conservative choice of power transformation.  In interactive mode, the user can look at the tightness of the scatter of
% points in the spread-vs-level plot and come to a sensible decision on the need to transform.
%
% fprev.  Fprev is included as an optional input argument to facilitate checking and modifying previously transformed series.
% The scenario is chronology development, in which an earlier step was to automatically transform all the ring widths. In 
% re-fitting, the user can leaf through the series one by one, evaluating suitability of each automatic transformation and 
% changing it if needed (see function powtran.m).  Re-fitting is indicated by existing data in fprev.  If fprev is not empty, 
% and if kopt(1)==4 (interactive fitting), the "initial fit" is the previous fit specified by fprev rather than the transformation 
% specified by the slope of the spread-vs-level plot.  If fprev is not empty and the fit mode is automatic, the previous 
% transformation information in fprev is ignored. 


close; % 

% Check fprev to find out whether this is a re-fitting or original fitting
% If a refit, unload old fit information
if isempty(fprev) | kopt(1)~=4;
    refit=0; % new fitting
else
    if ~iscell(fprev);
        error('fprev must be cell, if not empty');
    else;
        [mtemp,ntemp]=size(fprev);
        if mtemp~=1 | ntemp~=6;
            error('fprev must be 1 x 6');
        else; % unpack fprev
            cc=fprev{1};
            pp=fprev{2};
            aa=fprev{3};
            bb=fprev{4};
            eqnstrprev=fprev{5};
            khowprev=fprev{6};
            refit=1; % flag as a re-fit
        end;
    end;
end;


% Check X, and store as time series and year 
[mX,nX]=size(X);
if mX <30 | nX>2,
	error('X should be a 2-col Matrix, row size at least 30')
end
x=X(:,2);
yrx=X(:,1);
if any(diff(yrx)~=1);
    error('Year increment not 1 everywhere');
end;
if any(x<0);
    error('x must be non-negative');
end;
xmean=nanmean(x);
xsd=nanstd(x);

% x may possibly have NaNs.  Compute number of non-NaN observations
mX1 = sum(~isnan(x)); % sample size, ignoring any NaN

% Check segment length 
if lenseg>(mX1/3) | lenseg<10;
    error('lenseg must be at least 10, and less than 1/3 the number of non-NaN observations');
end;

% Check overlap
if mlap>=lenseg | mlap<0;
    error('mlap must be non-negative and less than lenseg');
end;

% Check options
[mkopt,nkopt]=size(kopt);
if mkopt~=1 | nkopt~=1;
    error('kopt must be 1 x 1');
end;
if kopt(1)<1 | kopt(1)>4; error('kopt(1) must be between 1 and 4');  end;

% Unpack text
id=txtin{1};
xlab1=txtin{2};

% Handle optional automatic log10 fit
if kopt(1)==1; % automatic log10
    khow=1;
    p=0;
    % Shift series if needed so that no zero values
    [ytemp,xinc]=shiftup1(x);
    clear ytemp; % unneeded
    % Transform
    c=xinc;
    [a,b,z]=subfn1(x,c,p);
    Y=[yrx z];
    eqnstr='y=log(x)';
    f={c,p,a,b,eqnstr,khow}; % put data in cells
    return;
end;


% SPREAD AND LOCATION (LEVEL) STATISTICS

% Shift series if any zero values
[x1,xinc]=shiftup1(x);

% Call pullseg2 to get a matrix whose columns are segments of the original time series. These segments can be
% overlapping.  The length of the segments is lenseg and the offset is mlap.  
[yry,Y]= pullseg2(x1,yrx,lenseg,mlap);
yry=yry'; % year vector corresponding to END years of the segments

% Compute the hinges for the segments. Hinges are min, lower quartile, median, upper quartile, max
H = hinge(Y); % each column has hinges for a segment

% Compute the approx interquartile range (iqr) and median for each segment
Q = (H(4,:) - H(2,:))'; % iqr
M = (H(3,:))'; % medians

% Log transform the segment medians and iqr's
u=log10(M);
v=log10(Q);
nsize=size(u,1); % number of segments

% Compute slope of log(iqr) vs log(median), and confidence interval for slope
U=[ones(nsize,1) u];
bslope= U\v;  % slope is bslope(2)
vhat=U * bslope; % fitted straight line
[B,BINT,R,RINT,STATS] = REGRESS(v,U,.01); % Get slope by regress function-- redundant, but want confid interval

% Compute power and its 95% confidence interval
p=1-bslope(2); % indicated power as function of slope
pint=1-BINT(2,:); % 95% interval around computed power, with value of 1.0 most conservative(no transformation
pstr=sprintf('%7.3f',p);

% Store initial computed-p info
pcalc=p;
pstrcalc=pstr;
pintcalc=pint;


% Store "initial fit" info
if refit==0; %  not a re-fitting
    p0=pcalc;
    p0int=pintcalc;
    p0str=sprintf('%7.3f',p0);
    eval(['eqnstr0=''y=x^{'  p0str '}'';']);
else; % This is a re-fitting
    p0=pp;
    p0str=sprintf('%6.3f',p0);
    eqnstr0=eqnstrprev;
end;


% AUTOMATIC FITTING USING SLOPE OF SPREAD VS LEVEL
if kopt(1)==2 | kopt(1)==3; 
    if kopt(1)==3; % power computed 1-slope
        khow=3;
        % p0, p0str, eqnstr0, p0int as already computed are applicable
        p=p0;
        pint=p0int;
        pst=p0str;
        eqnstr=eqnstr0;
        xinc=0;
        if p>0;
            % No shifting of zero values of x needed
        elseif p<0;
            % Need to add small amount to x if any values of x= 0
            [xtemp,xinc]=shiftup1(x);
            clear xtemp; % unneeded
        end;
        c=xinc;
        [a,b,z]=subfn1(x,c,p);
        eval(['eqnstr=''y=x^{' pstr '}'';']);
        f={c,p,a,b,eqnstr,khow}; % put data in cells
        Y=[yrx z];
        return;
        
    elseif kopt(1)==2; %  common power most compatible with slope of spread vs level
        khow=2;
        p=p0; % the power from slope of spread vs level
        pint=p0int;
        % Compute shift parameter that might be needed depending on the power
        [xtemp,xinc]=shiftup1(x);
        clear xtemp;
        c=xinc;
        % Set common powers 
        Pa=[3  2  1 2/3  1/2  1/3  0  -(1/3) -(1/2) -(2/3)  -1 -2]; % 12 common powers (0 corresp to log10, 1 to no transform)
        strPa={'y=x^{3}','y=x^{2}','y=x','y=x^{2/3}',  'y=x^{1/2}','y=x^{1/3}',...
                'y=log(x)','y-x^{-1/3}','y=-x^{-1/2}','y-x^{-2/3}','y=-x^{-1}','y=-x^{-2}'};
        if p>=1;
            Pb = [3 2 1]; % candidates
            strPb={'y=x^{3}','y=x^{2}','y=x'};
            PBo = [3 2 1]; % order of preference (1 most preferred)
            % Find common powers in 95% CI of computed power
            L1 = Pb>=min(pint) & Pb<=max(pint);
            if ~any(L1); % Not any common powers in 95% CI of computed power; select nearest common power
                L2=Pb<=p; % common powers more conservative than p
                if sum(L2)==1; 
                    ithis=find(L2);
                    
                else; % more than one candidate
                    ithis=find(L2);
                    ithis=ithis(1); % nearest
                end;
            else; % At least one common power in 95% CI of computed p
                if sum(L1)==1; % just one common power in the 95% CI
                    ithis=find(L1);
                else; % more than one common power in the 95% CI; Pick most conservative,even if not nearest
                    ithis = max(find(L1));
                end;
            end;
            p=Pb(ithis);
            eqnstr=strPb{ithis};
        elseif p<1;
            Pb=[ 1 2/3  1/2  1/3  0 -(1/3)   -(1/2)  -(2/3)   -1 -2];
            strPb={'y=x','y=x^{2/3}','y=x^{1/2}','y=x^{1/3}','y=log(x)','y=-x^{-1/3}','y=-x^{-1/2}','y=-x^{-2/3}','y=-x^{-1}','y=-x^{-2}'};
            Pbo=[1 2 3 4 5 6 7 8 9 10]; % order peference
            % Find common powers in 95% CI of computed power
            L1 = Pb>= min(pint) &  Pb<= max(pint);
            if ~any(L1); % Not any common powers in 95% CI of computed power; select nearest common power
                L2=Pb>=p; % common powers more conservative than p
                if sum(L2)==1; 
                    ithis=find(L2);
                else; % more than one candidate
                    ithis=find(L2);
                    ithis=max(ithis); % nearest
                end;
            else; % At least one common power in 95% CI of computed p
                if sum(L1)==1; % just one common power in the 95% CI; choose it
                    ithis=find(L1);
                else; % more than one common power in the 95% CI; Pick most conservative,even if not nearest
                    ithis = min(find(L1));
                end;
            end;
            p=Pb(ithis);
            eqnstr=strPb{ithis};
        elseif p==1; % no transform needed
            z=x; c=0; a=0; b=1; 
            eqnstr='y=x';
        end;
        [a,b,z]=subfn1(x,c,p); % transform & match
        % Transform
        f={c,p,a,b,eqnstr,khow}; % put data in cells
        Y=[yrx z];
        return;
    end
end;


%--  INITIAL SUBPLOTS

% Time series of original variable
subplot(2,2,1);
hp1=plot(yrx,x,[min(yrx) max(yrx)],[xmean xmean]);
ylabel(xlab1);
if isempty(id);
    id='Original Series';
end;
strtit1=['Time Series Plot: ' id];
legend('Series','Mean');
title(strtit1);

% Scatterplot of log(irq) vs log(median)
subplot(2,2,2);
hp2=plot(u,v,'*',u,vhat);
xlabel('Log(Median)');
ylabel('Log(iqr)');
set(gca,'XLim',[min(u) max(u)],'YLim',[min(v) max(v)]);
title('Spread vs Level for raw series, x');
pltext3(.05,.92,10,['Computed power, p = ' pstrcalc]);


% Time series plot of log(spread) and log(median)  
subplot(2,2,3);
[ax3,hp3left,hp3right]=plotyy(yry,u,yry,v);
xlabel(['End Year of ' num2str(lenseg) '-Year Period']);
set(hp3right,'LineStyle','--');
set(get(ax3(1),'YLabel'),'String','Log(Median)');
set(get(ax3(2),'YLabel'),'String','Log(iqr)');
legend('Median');
title('Time Changes in Spread and Level');
set(gca,'XGrid','on');

% Transform data according to p from initial transformation; intial depends on refit.  If refit==0, 
% initial is slope computed from spread vs level.  If refit==1, initial is previous transformation as 
% passed from calling function
p=p0;
if p>0;
    xinc=0;
elseif p<0;
    [xtemp,xinc]=shiftup1(x);
    clear xtemp;
end;
xinc0=xinc;
c=xinc;
[a,b,y]=subfn1(x,c,p);
y0=y;

% Shift y if needed to avoid inf values in logs
[yshift,yinc]=shiftup1(y);
clear yinc;

% Get the overlapping segments in Y;
[yrA,A]= pullseg2(yshift,yrx,lenseg,mlap);
yrA=yrA';
% yrA is the ending year of segments that are cols of A

% Compute the hinges for the segments
H = hinge(A);

% Compute the approx interquartile range for each segment
QA = (H(4,:) - H(2,:))';
MA = (H(3,:))'; % medians

% Log transform the segment medians and iqr's
u=log10(MA);
v=log10(QA);
nsize=size(u,1); % number of segments

% Compute slope
U=[ones(nsize,1) u];
bslope =U\v; 
% Get straight line fit
vhat=U * bslope;
% Compute power
p=1-bslope(2);


% Plot spread vs level, initial
subplot(2,2,4);
plot(u,v,'*',u,vhat);
xlabel('Log(Median)');
ylabel('Log(iqr)');
set(gca,'XLim',[min(u) max(u)],'YLim',[min(v) max(v)]);
p=p0;
pstr = sprintf('%7.3f',p);
if p<0;
    eval(['eqnstr=''y=-x^{' pstr '}'';']); 
elseif p>0;
    eval(['eqnstr=''y=x^{' pstr '}'';']); 
end;
title(['Spread vs Level for ' eqnstr0]);



%*********  Manual selection of transformation 

kwhile1=1;
kfirst=1;
while kwhile1;
    
    kmen2=menu('Choose one',...
        'Specify power',...
        'Log10 transform',...
        'No transformation',...
        'p=2/3',...
        'p=1/2',...
        'p=1/3',...
        'p=-1/3',...
        'p=-1/2',...
        'p=-2/3',...
        'p=-1',...
        'p=-2',...
        'Revert to initial transformation',...
        'Accept current transformation and exit while loop');
    
    if kmen2==1; % try new value for p
        khow=5;
        prompt={'Enter power, p:'};
        def={num2str(p)};
        dlgTitle='Specified value of power for tranformation';
        lineNo=1;
        answer=inputdlg(prompt,dlgTitle,lineNo,def);
        pstr=answer{1};
        p=str2num(pstr);
        if p0~=0; % if intial transform not log
            
            if abs((p-p0)/p0) <0.01;
                khow=4;
            else;
                khow=5;
            end;
        else; % initial transform was log; if current transform log, will not reach here
            
        end;
    
        x1=x;
        if p==0;
            close all; 
            error('p=0 is invalid power');
        end;
                
        if p<0; % negative power (inverse family)
            % zero to neg power is -inf; must shift input to avoid zeros
            [xtemp,xinc]=shiftup1(x);
            clear xtemp;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=-x^{' pstr '}'';']);
        else;
            xinc=0;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=x^{' pstr '}'';']);
            
        end
        kfirst=0;
    elseif kmen2==2; % log transform
        p=0;
        [xtemp,xinc]=shiftup1(x);
        clear xtemp;
        if p0==0;
            khow=4;
        else;
            khow=5;
        end;
        [a,b,y]=subfn1(x,xinc,p);
        eqnstr='y=log(x)';
        pstr=' ';
        kfirst=0;
    elseif kmen2==3; % No transformation
        if p0==1;
            khow=4;
        else;
            khow=5;
        end;
        y=x;
        p=1;
        pstr='1';
        eqnstr='y=x';
        kfirst=0;
        xinc=0;
    elseif kmen2==4; % p=2/3
        khow=5;
        p=2/3;
        pstr='2/3';
         
        x1=x;
        xinc=0;
        % Transform and match
        [a,b,y]=subfn1(x,xinc,p);
        eval(['eqnstr=''y=x^{' pstr '}'';']);
        kfirst=0;
    elseif kmen2==5; % p=1/2
        khow=5;
        p=1/2;
        pstr='1/2';
        
        x1=x;
        xinc=0;
        % Transform and match
        [a,b,y]=subfn1(x,xinc,p);
        eval(['eqnstr=''y=x^{' pstr '}'';']);
        kfirst=0;
    elseif kmen2==6; % p=1/3
        khow=5;
        p=1/3;
        pstr='1/3';
        
        x1=x;
        xinc=0;
        % Transform and match
        [a,b,y]=subfn1(x,xinc,p);
        eval(['eqnstr=''y=x^{' pstr '}'';']);
        kfirst=0;
    elseif kmen2==7; % p=-1/3
        khow=5;
        p=-1/3;
        pstr='-1/3';
        x1=x;
        if p<0; % negative power (inverse family)
            % zero to neg power is -inf; must shift input to avoid zeros
            [xtemp,xinc]=shiftup1(x);
            clear xtemp;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=-x^{' pstr '}'';']);
        else;
            xinc=0;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=x^{' pstr '}'';']);
            
        end;
        kfirst=0;
        
    elseif kmen2==8; % p=-1/2
        khow=5;
        p=-1/2;
        pstr='-1/2';
        x1=x;
        if p<0; % negative power (inverse family)
            % zero to neg power is -inf; must shift input to avoid zeros
            [xtemp,xinc]=shiftup1(x);
            clear xtemp;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=-x^{' pstr '}'';']);
        else;
            xinc=0;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=x^{' pstr '}'';']);
            
        end;
        kfirst=0;
    elseif kmen2==9; % p=-2/3
        khow=5;
        p=-2/3;
        pstr='-2/3';
        x1=x;
        if p<0; % negative power (inverse family)
            % zero to neg power is -inf; must shift input to avoid zeros
            [xtemp,xinc]=shiftup1(x);
            clear xtemp;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=-x^{' pstr '}'';']);
        else;
            xinc=0;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=x^{' pstr '}'';']);
            
        end;
        kfirst=0;
    elseif kmen2==10; % p=-1
        khow=5;
        p=-1;
        pstr='-1';
        x1=x;
        if p<0; % negative power (inverse family)
            % zero to neg power is -inf; must shift input to avoid zeros
            [xtemp,xinc]=shiftup1(x);
            clear xtemp;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=-x^{' pstr '}'';']);
        else;
            xinc=0;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=x^{' pstr '}'';']);
            
        end;
        kfirst=0;
    elseif kmen2==11; % p=-2
        khow=5;
        p=-2;
        pstr='-2';
        x1=x;
        if p<0; % negative power (inverse family)
            % zero to neg power is -inf; must shift input to avoid zeros
            [xtemp,xinc]=shiftup1(x);
            clear xtemp;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=-x^{' pstr '}'';']);
        else;
            xinc=0;
            % Transform and match
            [a,b,y]=subfn1(x,xinc,p);
            eval(['eqnstr=''y=x^{' pstr '}'';']);
            
        end;
        kfirst=0;
    elseif kmen2==12; % Revert to initial transformation
        khow=4;
        p=p0;
        pstr=p0str;
        eqnstr=eqnstr0;
        xinc=xinc0;
        y=y0;
        kfirst=0;
    elseif kmen2==13; %Accept current transformation and exit while looop
        if kfirst==1;
            khow=4;
            p=p0;
            xinc=xinc0;
            eqnstr=eqnstr0;
            [a,b,y]=subfn1(x,xinc,p);
            kfirst=0;
        end;
        kwhile1=0;
    end;
    
    % Shift to avoid log of zero
    [yshift,yinc]=shiftup1(y);
    
    % Get the overlapping segments in Y
    [yrA,A]= pullseg2(yshift,yrx,lenseg,mlap);
    yrA=yrA';
    % yrA is the ending year of segments that are cols of A
    
    % Compute the hinges for the segments
    H = hinge(A);
    
    % Compute the approx interquartile range for each segment
    QA = (H(4,:) - H(2,:))';
    MA = (H(3,:))'; % medians
    
    % Log transform the segment medians and iqr's
    u=log10(MA);
    v=log10(QA);
    nsize=size(u,1); % number of segments
    
    % Compute slope
    U=[ones(nsize,1) u];
    bslope= U\v; 
    vhat=U*bslope;
    
    subplot(2,2,4);
    plot(u,v,'*',u,vhat);
    xlabel('Log(Median)');
    ylabel('Log(iqr)');
    
    set(gca,'XLim',[min(u) max(u)],'YLim',[min(v) max(v)]);
    title(['Spread vs Level for  ' eqnstr]);
    
end % while kwh1


%-------  MATCHING

% Store transformed data and method flag
c=xinc;
[a,b,z]=subfn1(x,c,p);
f={c,p,a,b,eqnstr,khow}; % put data in cells
Y=[yrx z];
% END OF MAIN FUNCTION


function [a,b,z]=subfn1(x,c,p)
% subfn1: subfunction of sprdloc.m; computes parameters for matched transformation
% [a,b,z]=subfn1(x,c,p);
% Last revised 12-29-00
%
% Given a power time series x that had been transformed by a power or log transform to a new
% variable z, "matches" transformed variable z with the original series x such that the
% median and minimum are the same for the original and transformed series.  It is assumed that:
%  1) transform was power or log base 10, or none (p=1)
%  2) arbitrary constant c might have been added to x before transformation
% The matching is achieved by computing the parameters a and b of 
%  z = a + bT(x+c)  
% such that the minimim and median of z are equal to the minimum and median of x 

xmin=nanmin(x);
xmed=nanmedian(x);

if p==1;
    a=0; b=1; z=x;
    return;
elseif p==0;
    Tmin = log10(xmin+c);
    Tmed = log10(xmed+c);
    T=log10(x+c);
elseif p>0;
    Tmin = (xmin+c) .^ p;
    Tmed = (xmed+c) .^p;
    T = (x + c) .^ p;
else
    Tmin = -(xmin+c) .^ p;
    Tmed = -(xmed+c) .^ p;
    T = -(x + c) .^ p;
end;


% Compute parameter b
b = (xmed-xmin) / (Tmed - Tmin);

% Compute parameter a
a = xmed - b * Tmed;

z = a  + b * T;
Lneg=z<min(x);
if any(Lneg);
    z(Lneg)=min(x);
end;


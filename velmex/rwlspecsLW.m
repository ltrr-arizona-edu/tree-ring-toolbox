function rwlset=rwlspecsLW(rwlset,XEL);
% rwlspecsLW: interactive design of rwl specs for latewood width
% rwlset=rwlspecsLW(rwlset,XT,XEL);
% Last modified 2006-10-11
%
% Interactive design of rwl specs for latewood width.  
% Latewood width in conifers of NAM region typically lose signal for summer rainfall after
% first 150-200 years of growth of tree.  Trimming is therefore almost always necessary except
% for very young trees.
%
%*** INPUT 
%
% rwlset.  structure of .rwl file sets, with fields:
%   .name {j}s    short name or code of rwlset (e.g., padwt1}
%   .describe{j}   cell array of text description of the jth .rwl set
%   .trimall{j}    (1 x 2)i specified time coverage of data in the rwlset (all series truncated to this)
%   .idnames{j}  cell array of id names in jth set
%   .trimeach{j}  cell array of specified start and end years of spans for indiv series in the rwlset
% XEL structure with EW and LW data
%
%
%*** OUTPUT
%
% rwlset.  revised structure with an additional rwl set
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% .describe{}  rwlset.describe{1} holds the description for the first id set, or is [] if no
%   sets defined.  rwlset.describe{1} is itself a cell array of strings
% .idnames{}  rwlset.idnames{1} holds a cell array of series ids for series to be included in the
%   rwl set .  For example rwlset.idnames{1} ={'pdf04a','pdf05a'}
%
%

close all;
% Control over the window of "youth" considered by default to be useful to retain 
prompt={'Enter number of years:'};
name='By default, use only first n years of ringwidth series';
numlines=1;
defaultanswer={'200'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
ndef = str2num(answer{1});

p1 = 0.70 ; % will use spline with 50% response at 0.50 time sample length
 

id2 = XEL.id; 
num2=length(id2); % number of XEL ids
if num2==1 & isempty(id2{1});
    num2=0;
end;

% -- HOW MANY RWLSETS NOW EXIST

if isempty(rwlset.name{1});
    nsets=0;
else;
    nsets=length(rwlset.name);
end;
nextset=nsets+1; % next set, if creating new set
kmode='New'; % this function only generates a new rwlset

% -- WHAT TYPE OF SERIES WILL THE .RWL BE USED FOR -- THIS DETERMINES THE SET OF CANDIDATE SERIES
%  TO ICNLUDE

kmen =2; % the target is Earlywood-Latewood                   ');

ftype='ew';

% Build default menus
yesx = cellstr(repmat('-Y',num2,1));
nox = cellstr(repmat('-N',num2,1));
idset=id2;
numset=length(idset);
grpyx = cellstr([char(idset) char(yesx)]);
grpnx = cellstr([char(idset) char(nox)]);



%--- GET AND STORE FIRST AND LAST YEARS OF ALL LATEWOOD SERIES, AND BUILD STRING MENU TO SERIES

yrgo = repmat(NaN,num2,1); % to store start year of ring width series
yrsp = repmat(NaN,num2,1); % to store end year of ring width series
Lkeep = logical(ones(num2,1)); % flag to include the series in the rwlset
yron = repmat(NaN,num2,1); % to store start year of period of "good" data
yroff = repmat(NaN,num2,1); % to store end year of period of "good" data

T1=cell(1,num2); % to hold what will become rwlset{j}.trimeach;
T2 = cell(1,num2); % to hold the idnames

for n = 1:num2; % loop over LW series
    
    T1{n}=[];
    T2{n}=id2{n};
    
    X= XEL.data{n};
    yr = X(:,1);
    yrgo(n)=yr(1);
    yrsp(n)=yr(end);
    x = X(:,3);
    
    xdown = 0;
    xup = max(x)+0.05*(max(x)-min(x));
    
    % Start and end year of default youth window
    yr1 = yr(1);
    yr2 = yr(1)+ndef-1; 
    
    xpatch = [yr1 yr1 yr2 yr2 yr1];
    ypatch = [xdown xup xup xdown xdown];
    
    
    % Fit a cubic smoothing spline to the youth portion: g and yrg
     per = p1 * length(x);
     p = splinep(per,0.50);
     L1=yr>=yr1 & yr<=yr2;
     w1 = x(L1);
     yrw1 = yr(L1);
     L2 = isnan(w1); % mark any NaNs (if any, these are internal)
     w1temp = w1; % working copy of ring width time series
     w1temp(L2)=[]; % remove any NaNs from working copy of rw series
     yrw1temp=yrw1;
     yrw1temp(L2)=[]; % year vector for w1temp
     % Call csaps to fit the spline
     g = csaps(yrw1temp,w1temp,p,yrw1) ; % cv of spline growth curve, including values at any internal NaNs
     %   in the original ring widths w1
     
     yrg= yrw1;
     gc = max(g)/10;
     LL1 = g<gc;
     LL2 = g<0.10;
     LL3 = LL1 | LL2;
     
     i1 = find(LL3);
    
     if ~isempty(i1);
         yrc = yrg(i1(1));
         gc = g(i1(1));
     else;
         yrc = yr2;
         gc = g(yrg==yr2);
     end
     clear L1 w1 L2 yrw1 w1temp yrw1temp;
   
    figure(1);
    hp1=patch(xpatch,ypatch,[1 0.8 0.8]);
    set(gca,'XLim',[yr(1)-5  yr(end)+5]);
    hold on;
    
    hp2=plot(yr,x,yrg,g); 
    
    hline=line([yr(1) yr(end)],[0.05 0.05]);
    set(hline,'Color',[1 0 0],'Linewidth',2);
    
    
    grid on;
    title(id2{n});
    hold off;
    [cL,cB,cW,cH]=figsize(0.8,0.5);
    set(gcf,'Position',[cL cB cW cH]);
    
    kmen5=menu('Choose','Accept series with default ''Youth'' window',...
        'Accept, but click on end year for ''Youth Window''',...
        'Reject series');
    
    if kmen5==1;
        yron(n)=yr1;
        yroff(n)=yr2;
    elseif kmen5==2;
        yron(n)=yr1;
        [xpt,ypt]=ginput(1);
        yroff(n)=floor(xpt); 
    else
        yron(n)=NaN;
        yroff(n)=NaN;
        Lkeep(n)=0;
    end
    
    close(1);
    T1{n}=[yron(n) yroff(n)];
    
    
end; % % loop over LW series

%-- Name the new rwlset
prompt={'Enter short name for rwl set','Enter description:'};
def={'DefaultSet','Default set of series to send in .rwl file to ITRDB'};
dlgTitle='Name and description for this rwlset';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
rwlset.name{nextset}=answer{1};
rwlset.describe{nextset}=answer{2};

% Store the creation date
rwlset.when{nextset}=date;

% specify no restriction on the time coverage of the rwlset
rwlset.trimall{nextset}=[];

rwlset.idnames{nextset}=T2(Lkeep); % The ids of cores in this set

rwlset.trimeach{nextset}=T1(Lkeep); % The trim period of ring widths for the rwl set

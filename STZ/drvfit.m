function crvfit
% EDITING TO:
%  * prompt for diff vs ratio fitting; if difference, need to check that A exist--DONE
%  * revise load and save to take into account tranform.m contributions. MAY have  A,tstatus,Fwhen
%  * revise to build and handle fitting history a la Fwhen
%  * Somehow record whether index from diff or ratio.  --   DONE -- AS S(:,6)
%  * If diff,
%       compute & temp store transformed analog of X; operate on that
%       est age trend on X by BOTH diff and ratio methods
%       Scale diff index so that mean 1.0 and variance eq to that of ratio index
%       Store scaling info, G, etc
%
% 
% crvfit: interactive fitting of curves to detrend ring-width series
% crvfit
% Last revised 1-8-01
%
% Select curve fit for detrending ring width.  Fit the detrending curve. 
% Store the information on the curve choices and the fitted trend lines.
%
%*** IN ****************************************
%
% User prompted for name of .mat file storing the ringwidth data
% and associated years and names.  Ringwidth data previously put in
% this file with rwlinp.m
%
%*** OUT ***********************************
%  

% User prompted for name of .mat file storing curve-fit information. Can
% store in a new file, or add to the .mat file that holds the input
% This file may will contain new or revised copies of 
%  nms, X, yrs
%  S 
%  G
%  Fwhen {8 x 4} history of chronology development (date, hour, minute that various programs rin
%   {1,:} rwlinp
%   {2,:} tranform
%   {3,:} crvfit
%   {4,:} corei
%   {5,:} treei
%   {6,:} sitei
%   The four cells for each row hold
%       function : e.g., "rwlinp"
%       date&time: e.g., Jan 01 01  10 3 == Jan 1 2001, hour 10, minute 3
%       input file:  e.g., "ste.rwl"
%       output file:  e.g., "ste.mat"
%   {7,:} itrdbcrn
%   {8,:} crntbl
%  
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED (FROM C:\MLB\STZ UNLESS NOTED)
%
% blocdat
% cfnegx
% erchk
% monotspl
% cfstrln1
% cfspl
% negxpk1 (c:\mlb)
% splinep
% cfmean
%
%*** TOOLBOXES NEEDED
%
% STATS
% OPTIM
% SPLINES
%
%*** NOTES
%
% The storage matrix S holds the curve fit information, one line per core
%
% col 1 -- sequence number of core, or zero if no curve yet fit
% col 2 -- curve fit option for "first detrending"
%   1 modified negative exponential
%   2 straight line
%   4 horizontal line at mean
%   9 spline
% cols 3-5  parameters for curve fits; content varies according to type of curve
%   1 (neg exp): 3:5 are k, a, b
%   2 (str line):  3 is a, 4 is b, in eqn y=a*x + b;  5 is zero
%   4 (horiz line at mean):  3 is mean, 4 and 5 are zero
%   9 (cubic spline):  3:5 are p, per, amp, where p is the spline parameter,
%		per is the period in years and amp is the amplitude of freq response 
%		at period per
% col 6:  Ratio vs difference for index computation
%       ==1 ratio
%       ==2 difference
% col 7: scale parameter to adjust difference-method indices so that the index mean is 1.0 and 
%   the index standard deviation is equal to the standard deviation of a ratio-index derived using the same 
%   fitted growth curve as the difference index (see notes)
%
% cols 6-9 reserved for future use.  Previously these slots used to store
%  "second detrending" data analagous to cols 2-5.
% col 10, 11 -- start, end year of period for fitting curve
% col 12 -- starting row index in sov X for year in col 10
%
% User defines endpoint of segment to be fit
% Detrending curve is fit to the segment
% User can block out sub-segments and re-fit smoothed curve to
% 	remaining data
% User accepts or revises detrending
%
%
%*** NOTES
%
% Scale parameter of difference-method index.  Say the growth curve is g and the matched, tranformed ring width 
% is w.  The difference index is then D=w-g and the ratio index is R=w/g. Let d=D-Dmean and r=R-Rmean be the mean-subtracted 
% difference and ratio indices, and std(D) and std(R) be the standard deviations of those indices.  The shifted, scaled 
% difference index is computed as D1 = 1.0 + [std(D)/std(R)]*d = b * d.  The value for b is stored as S(:,7).  If index 
% method is ratio, S(:,7) is set to zero.


% Close any open windows
close all

% Hardcode strings for 9 curve types
strfit={'Mod. Negative Exponential','Least Squares Straight Line',' ','Horizontal Line at Mean',...
        ' ', ' ', ' ', ' ', 'Spline'};

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

% Fit history
Lwhen=exist('Fwhen','var')==1;
if ~Lwhen;
    Fwhen = cell(6,4);
end;
Fwhen{3,1}='crvfit'; % function
Fwhen{3,3}=flmat; % infile

clear flmat


%***********************************************************************

%--- Prompt for ratio vs difference detrending
kmen2=menu('Choose mode for detrending',...
    'Ratio',...
    'Difference');
switch kmen2;
case 1;
    detmode=1; % ratio
    strdet='Ratio';
case 2; %detmode='Difference';
    detmode=2; % difference
    strdet='Difference';
    if exist('A','var')~=1;
        error('Difference method, but no A with stored transformation parameters');
    else; % check that ALL cores had been transformed
        if any(tstatus==0);
            error('Not all core ring widths transformed yet');
        end;
    end;
end; % switch
  

%-- Prompt for all-the-same curve fits
autocrank=questdlg('Same curve-fit for all cores?');
if strcmp(autocrank,'Yes');
   kmen1=menu('Choose automatic curve fit for all series',...
      'Nonincreasing spline',...
      '0.8N spline',...
      'Horizontal Line',...
      'Negative-Slope Combo');
   if kmen1==1;
      curvetype='sni';
   elseif kmen1==2; 
      curvetype='s8';
   elseif kmen1==3;
      curvetype='HL';
  elseif kmen1==4;
      curvetype='nsc'; % Negative-slope combo
      % non-increasing spline first try
      % if slope positive, pick horiz line at mean
      % if if ni spline OK but %N > 2N, take SL
   end;
   
   

   
   
else;
    
    %  Prompt for specs of robust-mean low frequency index
   prompt={'Enter the wavelength at which amplitude of response 0.50:'};
   def={'100'};
   dlgTitle='0.50 Response Wavelength of Spline for Robust-Mean Low Frequency Index';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
   plfi=str2num(answer{1}); % period of robust low-frequenc 
   
   %  Prompt for specs of robust-mean low frequency index
   kmen2=menu('Robust-mean low frequency index will be computed from (choose one):',...
       '100% of series overlapping test series',...
       'User specified precentage of overlap series',...
       'Abort');
   if kmen2==1; 
       pctuse=100;
       close all;
       fclose all;
       disp(['0.50 period of LFI = ' num2str(plfi) ' yr']);
       disp(['pctuse = ' num2str(pctuse) ]);
       
       
   elseif kmen2==2;;
       prompt={'Enter the percentage:'};
       def={'75'};
       dlgTitle='Percentage if overalapping series to use in computing robust mean low-frequency index';
       lineNo=1;
       answer=inputdlg(prompt,dlgTitle,lineNo,def);
       pctuse=str2num(answer{1}); % percentage of series to use.  The LONGEST overlapping series are given
       % precede nce. 
       close all;
       fclose all;
       isp(['0.50 period of LFI = ' num2str(plfi) ' yr']);
       disp(['pctuse = ' num2str(pctuse)]);
       
   else;
   end;
   
   
   %-- ROBUST LOW-FREQUENCY COMPONENT
   %
   % Will store the tsm of low frequency components for each series in H, with year vector yrH
   
   % Compute low frequency components for individual series
   
   
   % Must have fitted growth curves beforehand
   if ~exist('G','var')==1;
       error('No growth curves yet; go back and run in autocrank mode');
   end;
      
   % If have masked out cores, will not use those for the robust LFI
   if exist('cmask','var')==1;
       coremask=cmask;
   else;
       coremask = logical(ones(size(nms,1),1));
   end;
   
   % If called for difference indexing, will need parameter matrix A
   if detmode==1; % ratio detrending
       AA=[];
   else;
       if exist('A','var')~=1 ;
           error('Wanted difference index, but no parameter matrix A  created yet');
       else;
           AA=A;
       end;
   end;
   % Call subfunction to compute the individual series LFIs. Recall that coremask is zero to 
   % mask out series, that AA deals with difference indexing, and that plfi is the period 
   % at which the amplitude frequency response of the spline used to comput the LFI is 0.5
   [H,yrH,rath]=subfn01(X,yrs,S,G,AA,coremask,plfi);
      
   curvetype='null';
   kmen1=NaN;
end;

%-- Prompt for whether to use "blocking out".
if strcmp(autocrank,'No');
   blockuse=questdlg('Will you be blocking out segments in curve fitting?');
   blockuse=upper(blockuse); % change blockuse to upper case
   switch blockuse;
   case {'YES','CANCEL'};
       error('Blocking not yet checked out -- answer No to query on using blocking');
   end;
       
else; % autocrank
   blockuse='NO';
end;



% Get the number of cores, which equals the row size of  matrix
% of core ids
ns=size(nms,1);


% If no growth curve info yet in the storage file, initialize
% G
mX=size(X,1);  
if exist('G')~=1;
   G = repmat(NaN,mX,1);
end;


% Get information on whether particular cores have already been fit
% in a previous run of crvfit.m.  If this is a first run, initialize S and the
% string matrix fv. If S exists, col 1 of S will either be a sequence
% number or 0, depending on whether the core has been fit or not.
fvv='-n'; %  not fit yet
fvf='-f'; %  already fit 

% If S exists, get the info on whether cores have been fit and store in fv,
% If S does not exist, initialize as zeros
fv = repmat(fvv,ns,1); % initialize string matrix as if no cores yet fit
if exist('S'),
  LS1=S(:,1)~=0;  % non-zero in first col of S:  core had been fit
  numf = sum(LS1);  % number of cores already fit
  fv (LS1,:) = repmat(fvf,numf,1);
else; % If this is first run of crvfit.m, intialize matrix S as zeros
  S=zeros(ns,22);
end

% Make a string matrix with the following for each row:
%  * the sequential number of the core
%  * '-'
%  So, might get
%	1-
% 	2-
%  etc
nsmt = [int2str((1:ns)')  repmat('-',ns,1)];


ksw1=1; % while switch for continue working on this core

kfx=0.2; % x and y plotting positions for tsp of ringwidth and fitted curve
kfy=0.15;


%******************  WHILE OVER CORES  ********************

iauto=0;  % initialize sequence number for all-the-same fits

while ksw1==1;  % while working on this core
   iauto=iauto+1;
   nsms=[nsmt nms fv]; % string matrix with sequence number,
   % nsmt is sequential number plus -
   % nms is core id
   % fv is fit status: -f == "fit", '-n'="not yet fit" 
   
   % Select a core. 
   if strcmp(autocrank,'No');
      scid=menu('Core ID # ?',(cellstr(nsms))');  % sequence number of selected core
   else;
      scid=iauto;
   end;
   
   
   % Get ring-width series
   xv=X(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));
   
   % If difference method, re-compute transformed, matched ring width
   if detmode==2; % diff index
       cshift = A{scid,1}; % inc to add to rw before power tran
       ppower= A{scid,2}; % power of transformation
       acoef = A{scid,3}; % shift coef for matching
       bcoef = A{scid,4}; % mult coef for matching
       
       % Power tran
       if ppower==1; % moot power; no transformation
           xvtran=xv;
       elseif ppower==0; % log transform
           xvtran=log10(xv+cshift);
       elseif ppower>0; % positive power
           xvtran= (xv+cshift) .^ ppower;
       elseif ppower<0; % neg power
           xvtran= -(xv+cshift) .^ ppower;
       end;
       
       % Matching
       xvmatch=acoef + bcoef * xvtran;
       xv=xvmatch;
   else;
   end;
   
       
   % Get stored previously fitted growth curve -- or NaN if not yet fit
   gv=G(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1)); 
   
   
   % Build title string for curvfit
   % S(:,2) is the curve type: 1=NE, 2=   , 4=HL,  9=spline
   % For Splines: S(:,4) is length(yr) and 
   if exist('S','var')==1 & S(scid,2)~=0;
       
       strfita = strfit{S(scid,2)};
       if S(scid,2)==9 ; % if spline
           strfit1b = num2str(S(scid,4),'%4.0f');
           strfita = [strfita ', ' strfit1b ' yr ('];
           % Length of "fit" period
           nyrfit = S(scid,11)-S(scid,10)+1;
           fitratio = S(scid,4)/nyrfit;
           strfitc =[ num2str(fitratio,'%5.2f') 'N)'];
           strfita=[strfita strfitc];
           clear strfitb strfitc;
       end;
   else;
       strfita=' ';
   end;

   
   yrv=(yrs(scid,1):yrs(scid,2))';  % col vector of years for xv and gv
   
   if strcmp(autocrank,'No');
       
       % Pull the robust LFI computed on all cores excluding this
       Lr1 = yrH >= min(yrv) & yrH<=max(yrv);
       Lc1 = coremask;
       Lc2 = logical(ones(ns,1));
       Lc2(scid)=0;
       hkey = H(Lr1,scid);
       H1 = H(Lr1,Lc1&Lc2);
       % Further require that any series for the robust LFI start earlier than the test series
       Lc3=~isnan(H1(1,:));
       H1 = H1(:,Lc3);
       if size(H1,2)>1;
           h1 =   (nanmedian(H1'))';
       elseif size(H1,2)==1;
           h1=H1;
       else;
           h1=repmat(NaN,length(yrv),1);
       end;
       
       % Compute number of series in robust LFI each year
       if size(H1,2)>1;
           Nrobust =   (nansum(~isnan(H1')))'; % cv of number of series
       elseif size(H1,2)==1;
           Nrobust = ~isnan(H1);
       else;
           Nrobust=zeros(length(yrv),1);
       end;
       Ldecade=rem(yrv,20)==0;
       yrvrobust=yrv(Ldecade);
       Nrobust=Nrobust(Ldecade);
       strrobust=num2str(Nrobust);
              
       % Scale the robust LFI to same mean and variance of fitted growth curve
       if ~all(isnan(h1));
           %h1 = (h1 - mean(h1))/std(h1)  * (std(xv)* rath(scid)) + mean(h1);
       end;
       
   end;
   
   
   % Plot ring width or transformed, matched ring width
   if strcmp(autocrank,'No');
       hf0=figure('Units','normal','Position',[kfx,kfy,0.75,0.75]);
       if strcmp(fv(scid,:),'-n');
           plot(yrv,xv,'b');title(['Ring Width ',nms(scid,:)]);
           xlabel('Year');ylabel('Ring Width, (0.01 mm)');
       else;
           if all(isnan(h1));
               plot(yrv,xv,'b',yrv,gv,'m');
               title(['Ring Width ',nms(scid,:) '& curve fit: ' strfita]);
           else;
               haxis1=axes('Position',[0.13 0.33  0.80 0.63]); 
               hpplot1=plot(yrv,xv,'b',yrv,gv,'m');
               ylabel('Ring Width, (0.01 mm)');
               title(['Ring Width ',nms(scid,:) '& curve fit: ' strfita]);
               haxis2=axes('Position',[0.13 0.10 0.80 0.23]);
               hpplot2=plot(yrv,h1,yrv,hkey,[min(yrv) max(yrv)],[1 1]);
               xlabel('Year');
               ylabel('Index');
               set(hpplot2(1),'Color',[.7 .7 .7],'LineWidth',2);
               text(yrvrobust,h1(Ldecade),strrobust);
           end;
       end;
       
   end;
   
   
   xv1=xv; % store full-length ring-width series
   gv1=gv; % store full-length grwoth trend
   
   % If growth curve previously fit, allow to accept the old fit and move on
   if strcmp(autocrank,'No');
       if strcmp(fv(scid,:),'-f');
           keeper=questdlg('Accept previous fit');
       end;
       keeper=upper(keeper);
       if strcmp(keeper,'YES');
           ksw2=0;
       elseif strcmp(keeper,'NO');
           ksw2=1;
           gv=repmat(NaN,length(gv),1);
       else; % bail out and allow exportfig
           fclose('all');
           return;
           
       end;
   else; % autocrank mode
      ksw2=1;
      gv=repmat(NaN,length(gv),1);
   end;
   
      
         
   %**************  FIRST DETRENDING **************
   
   while ksw2==1;  % working on first detrending
      
      % Prompt for interval of series to fit detrending curve to
      if strcmp(autocrank,'Yes');
         nseg=4; % in autocrank mode, curve is fit to full length
      else;
         nseg=menu('Ends ?','Graphical input','Specify in #','Full Length');
      end;
      if nseg==1, % graphically point to ends of fit period
         jh1=msgbox('Click at two corner points of fit period',' ');
         pause(1);
         close(jh1);
         figure(hf0);
         [segx,segy]=ginput(2);
         [erflg,xind1,xind2]=erchk(yrv,min(segx),max(segx));
      elseif nseg==2, % specify start and end years of fit period
         prompt={'Enter start year: ','Enter end year: '};
         titledlg='Period to fit curve to';
         def={int2str(yrv(1)),int2str(yrv(length(yrv)))};
         LineNo=1;
         answer=inputdlg(prompt,titledlg,LineNo,def);
         yeargo=str2num(answer{1});
         yearstop=str2num(answer{2});
         % compute row indices in xv and yrv corresponding to the selected period
         [erflg,xind1,xind2]=erchk(yrv,yeargo,yearstop);
      else; % fit detrending curve to entire ring-width series
         xind1=1;
         xind2=length(yrv);
         erflg=0;
      end
      
      % Make sure specified year range for fit period valid
      if erflg==-1,
         close all
         fclose all
         error('Years out of range');
      end
      
      yrvn=yrv(xind1:xind2); % year vector for selected fit interval
      xvn=xv(xind1:xind2);   % ring-width data for the selected fit interval
      
      S(scid,10)=yrv(xind1); % Store start year for fit interval
      S(scid,11)=yrv(xind2); % Store end year for interval
      S(scid,12)=yrs(scid,3)+xind1-1; % store starting index for interval
      
      % Plot time series of ring width for selected fit interval
      if strcmp(autocrank,'No');
         hf1=figure('Units','normal','Position',[kfx,kfy,0.75,0.75]);
         plot(yrvn,xvn,'b');title(['RingWidth ',nms(scid,:)]);
         xlabel('Year');ylabel('RW, (0.01 mm)');
      end;
            
      param=repmat(NaN,1,8);  % Initialize parameters for curve fit
      
      % Prompt for curve fit options and Curve fitting
      % eventual curve options ['1';'2';'3';'4';'5';'6';'7';'8';'9'];
      nopt1=['1 - Neg Ex';'2 - SL    ';'4 - Mean  ';'9 - Spline'];
      
      % Set flag for type of detrending curve
      if strcmp(autocrank,'No');
         nfit = menu('Curve-fit option',(cellstr(nopt1))');
      else;
         if strcmp(curvetype,'s8') | strcmp(curvetype,'sni') | strcmp(curvetype,'nsc');
            nfit=4;
         elseif strcmp(curvetype,'HL');
            nfit=3;
         end;
      end;
      
      if nfit==1; % neg exp
         nfit1=1;
     elseif nfit==2; % straight line, any slope
         nfit1=2;
     elseif nfit==3; % horiz thru mean
         nfit1=4;
     else
         nfit1=9; % spline
     end
      strfita=strfit{nfit1};
      
      % Store curve-type for detrending
      S(scid,2)=nfit1;
      
      if nfit1==9; % Cubic smoothing spline
         if strcmp(autocrank,'No');
            splp=menu('p-option ?','%N and .5 Amp','#yrs and .5 Amp','%N and x Amp',...
               '#yrs and x Amp','Specify p','Non-Increasing');
         else; % autocrank mode
            if strcmp(curvetype,'s8');
               splp=1;
            elseif strcmp(curvetype,'sni') | strcmp(curvetype,'nsc');
               splp=6;
            end;
         end;
         
         if splp>7;
            jh1=jdisp('Please return to command window for input');
            pause(1);
            close(jh1);
         end;
                
         
         if splp==1; % you specify the wavelength with 0.5 AFR as fraction of series length
            if strcmp(autocrank,'No');
               prompt={'Enter the decimal fraction of series length:'};
               def={'.70'};
               dlgTitle='0.50 Frequency Response Wavelength';
               lineNo=1;
               answer=inputdlg(prompt,dlgTitle,lineNo,def);
               pper = str2num(answer{1});
               %pper=input('Period decimal fraction of the total series length = ');
            else; % autocrank
               pper=0.8;  % spline with 0.5 amp freq response at 0.8 sample length
            end;
            per=pper*length(xv);
            amp=0.5;
            p=splinep(per,amp);
         elseif splp==2; % For 0.5 AFR, specify wavelength in number of years
            prompt={'Enter the Number of Years:'};
            def={'128'};
            dlgTitle='0.50 Frequency Response Wavelength';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            per = str2num(answer{1});
            %per=input('Length of the period (Years) = ');
            amp=0.5;
            p=splinep(per,amp);
         elseif splp==3; % Specify AFR specified AFR at specified wavelength
            prompt={'Enter Desired wavelength as fraction of N:','Enter the AFR:'};
            def={'0.7','0.5'};
            dlgTitle='Desired \lambda as %N and AFR';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            pper = str2num(answer{1});
            amp = str2num(answer{2});
            
            %pper=input('Period as % of the total series length = ');
            per=pper*length(xv);
            %amp=input('Please enter the value of Amplitude = ');
            p=splinep(per,amp);
         elseif splp==4; % specify wavelength as number of years and AFR
            prompt={'Enter Desired wavelength (years):','Enter the AFR:'};
            def={'128','0.5'};
            dlgTitle='Desired \lambda in years, and AFR';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            per = str2num(answer{1});
            amp = str2num(answer{2});
           
            %per=input('Length of the period (Years) = ');
            %amp=input('Please enter the value of Amplitude = ');
            p=splinep(per,amp);
         elseif splp == 5;
            prompt={'Enter spline parameter:'};
            def={'1E-5'};
            dlgTitle='Spline Parameter';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            p=str2num(answer{1});
	         %p=input('Spline parameter, p = ');
            per=NaN;
            amp=NaN;
        elseif splp ==6; 
            if strcmp(autocrank,'No');
               [p,per]= monotspl(yrv,yrvn,xvn,length(xv),1);
            else;
               [p,per]= monotspl(yrv,yrvn,xvn,length(xv),2);
           end;           
           amp = 0.5;
       end;
         if strcmp(autocrank,'No'); % not in autocrank
             param(2:4)=[p per amp];
             cvx = (cfspl(p,yrv,yrvn,xvn))'; % Compute spline
             tstr1=['  SPL, p = ',num2str(p)];
         else; % autocrank mode
             cvx = (cfspl(p,yrv,yrvn,xvn))'; % Initially compute age curve by spline
             % If first value of age curve lower than last, have monotonic increase: call for HL
             if cvx(length(cvx))>cvx(1);
                 cvx = cfmean(yrv,yrvn,xvn); % Horizontal line through mean
                 tstr1='  HMN instead of +slope spline ';
                 param(2)=mean(xvn(~isnan(xvn))); % sample mean of the valid ringwidths
                 param(3:4)=NaN;
                 nfit=3;
                 nfit1=4;
                 S(scid,2)=nfit1;
             else;  % spline is nonincreasing; but period of 0.5 AFR > 2N.  Fit SL
                 if per > 2 * length(cvx);
                     hwarn = 0;  % unneeded for this call
                     [cvx,a,b] = cfstrln1(yrv,yrvn,xvn,hwarn);
                     param(2:3)=[a b];
                     param(4)=NaN;
                     tstr1 = '  SL neg slope ';
                     nfit=2;
                     nfit1=2;
                     S(scid,2)=nfit1;
                 else;  % accept the auto-fitted nonincreasing spline
                     param(2:4)=[p per amp];
                     tstr1=['  SPL, p = ',num2str(p)];
                 end;
             end;
         end;
      
      elseif nfit1==1; % Neg exp
         hwarn = 1;  % want to display warning dialog if neg exp wrong type
         [cvx,k,a,b] = cfnegx(yrv,yrvn,xvn,hwarn);
         param(2:4)=[k a b];
         tstr1 = '  NEG EXP ';
         
         % If you want to generate the curve in an outside program.
          % g(t) = k + a * exp(-b*t),
          %   where t is the shifted time variable t = yrvn-yrvn(1)+1
          %   In other words, t is same length as yrvn after 
          %   dropping NaNs
          
         
      elseif nfit1==2; % straight line, any slope
         hwarn = 0;  % unneeded for this call
         [cvx,a,b] = cfstrln1(yrv,yrvn,xvn,hwarn);
         param(2:3)=[a b];
         tstr1 = '  SL any slope ';
   
      elseif nfit1==4; % Horizontal line through mean
         cvx = cfmean(yrv,yrvn,xvn); % Horizontal line through mean
         tstr1='  HMN ';
         param(2)=mean(xvn(~isnan(xvn))); % sample mean of the valid ringwidths
      end
      
      S(scid,3:5) = param(2:4); % Store parameters for curve-fit
      
      % Make string for curve fit
      strfitb = strfit{S(scid,2)};
      if S(scid,2)==9 ; % if spline
          strfitc = num2str(S(scid,4),'%4.0f');
          strfitb = [strfitb ', ' strfitc ' yr ('];
          % Length of "fit" period
          nyrfit = S(scid,11)-S(scid,10)+1;
          fitratio = S(scid,4)/nyrfit;
          strfitc =[ num2str(fitratio,'%5.2f') 'N)'];
          strfitb=[strfitb strfitc];
          clear strfitc nyrfit fitratio ;
          
      end;
      
      % Plot trend-line superposed on ring width
      if strcmp(autocrank,'No');
         figure(hf1);hold off;
         plot(yrvn,xvn,'b',yrvn,cvx(xind1:xind2),'k');
         title(['Ring Width ',nms(scid,:), ' & Trend by ' strfitb ]);
         %title(['ERE ' strfita]);
         xlabel('Year');ylabel('Ring Width, (0.01 mm)');
      end;
      
      % Initialize blocking settings
      ksw3=1; % while control for blocking segments from use in curve fit
      nblk=1;
      sgc=1;
      
      %********** WHILE BLOCKING OUT INTERVALS ***********
      
      while ksw3==1 & nblk~=6 & strcmp(blockuse,'YES');
         % Blocking out specified data segments
         if sgc==1,
            xv1=xv;
         end
         [S,yrvn,xvn,xv1,nblk] = blocdat(hf1,scid,xind1,yrvn,...
            xvn,yrv,xv,xv1,sgc,S);
         if nblk==6, break; end
         figure(hf1);hold off;
         plot(yrvn,xvn,'b'); title(['Ring Width ',nms(scid,:)]);
         xlabel('Year');ylabel('Ring Width, (0.01 mm)');
         pause(1);
         ksw3=menu('Select','Block new segments ?','Continue');
         if ksw3~=1,
            % Curve fitting for the blocked data
            if nfit1==9; % Spline
               cvx = (cfspl(p,yrv,yrvn,xvn))';
            elseif nfit1==1; % Neg exp
               hwarn = 1;  % want to display warning dialog if neg exp wrong type
               [cvx,k,a,b] = cfnegx(yrv,yrvn,xvn,hwarn);
               param(2:4)=[k a b];
               S(scid,3:5) = param(2:4); % Store parameters
               tstr1 = '  NEG EXP ';
            elseif nfit1==2; % straight line, any slope
               hwarn=0;
               [cvx,a,b]=cfstrln1(yrv,yrvn,xvn,hwarn);
               param(2:3)=[a b];
               S(scid,3:4)=param(2:3);

            elseif nfit1==4; % Horizontal line through mean
               cvx = cfmean(yrv,yrvn,xvn);
               param(2)=nanmean(xvn); % sample mean of the valid ringwidths
               S(scid,3)=param(2);

            end
            
            strfitb = strfit{S(scid,2)};
            if S(scid,2)==9 ; % if spline
                strfitc = num2str(S(scid,4),'%4.0f');
                strfitb = [strfitb ', ' strfitc ' yr ('];
                % Length of "fit" period
                nyrfit = S(scid,11)-S(scid,10)+1;
                fitratio = S(scid,4)/nyrfit;
                strfitc =[ num2str(fitratio,'%5.2f') 'N)'];
                strfitb=[strfitb strfitc];
                clear strfitc nyrfit fitratio ;
                
            end;
            
            
            % Plot trend line superposed on ring width
            figure(hf1);hold off;
            plot(yrvn,xvn,'b',yrv(xind1:xind2),cvx(xind1:xind2),'k');
            title(['Ring Width ',nms(scid,:),tstr1 ' ' strfitb]);
            xlabel('Year');ylabel('Ring Width, (0.01 mm)');
         end
         sgc=sgc+1;
      end 			% End of ksw3 (blocking) while loop
      %***********************************************************
      
      
      % Allow user to change options on  dtrending
      if strcmp(autocrank,'No');
         ksw2=menu('Select One','Change curve-fit ?','Accept curve-fit'); 
         if ksw2==2,
            nq=logical(1);
         elseif ksw2==1;
            nq=logical(0);
         else
            error('Impossible value for ksw2');
         end
      else; % autocrank
         ksw2=2;
         nq=logical(1);
      end;
      
      if nq; % you accepted curve-fit
          
          % Ratio version of index
          rwin1ratio=xv1(xind1:xind2)./cvx(xind1:xind2);
          mnxratio=nanmean(rwin1ratio); % mean of ratio index
          rwin2ratio = rwin1ratio - mnxratio +1.0; % Adjusted to a mean of 1
          stdratio=nanstd(rwin2ratio); % std dev of ratio index
          
          % Difference version of index
          rwin1diff=xv1(xind1:xind2) - cvx(xind1:xind2);
          mnxdiff=nanmean(rwin1diff); % mean of diff version of index
          depdiff = rwin1diff-mnxdiff; % diff index as departures from mean
          stddiff = nanstd(rwin1diff); % std dev of diff index, before scaling
          sfact=(stdratio/stddiff);
          depdiff=depdiff * sfact; % scale departures so that std same as for ratio index
          rwin2diff = depdiff + 1.0; % shift to mean of 1
          
          % Store the appropriate (ratio vs diff ) index
          switch detmode;
          case 1; % ratio
              rwin2=rwin2ratio;
          case 2; % difference
              rwin2=rwin2diff;
              S(scid,7)=sfact; % scale factor to be multuplied by departures of diff index from mean before adding 1.0
              % to get final scaled diff index.
          end;
          
          
          % compute index as ratio of ring-width to fitted curve
         rwin1=xv1(xind1:xind2)./cvx(xind1:xind2);
         % Comnpute  mean of ratio ring-width to fitted curve 
         mnx=nanmean(rwin1);
         % compute index adjusted to mean 1.0
         rwin2 = rwin1 - mnx +1.0;
         
         % Plot index
         if strcmp(autocrank,'No');
            hf2=figure('Units','normal','Position',[kfx,kfy,0.75,0.75]);
            plot(yrv(xind1:xind2),rwin2,'b',[yrv(xind1) yrv(xind2)],[1.0 1.0],'r');
            strtemp=[' Tree-Ring Index by ' strdet ' Method'];
            title([nms(scid,:), strtemp]);
            
            xlabel('Year');ylabel('Index');
         end;
         
      else; % want to re-do curve fit
         ksw2=1;
      end
      
      % Store the growth trend
      gv(xind1:xind2)=cvx(xind1:xind2);
      G(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1))=gv; % Store fitted trend line
      
      % Update the Robust Low Frequency Index matrix
      H(:,scid)=NaN;  % sub NaN for existing stored low-frequency component of test series
      % recall that S(scid,10) has start year for fit period and S(scid,11) has end year for fit period
      % recall that rwin2 is the current version of the core index
      % Recall that yrH is the year vector for the robust low-freq index matrix H
      i1H = find(yrH==S(scid,10)); % index to row of H to store first value in
      i2H = find(yrH==S(scid,11)); % index ... last value
      yrthish = ((S(scid,10)):(S(scid,11)))';
      % Compute low-frequency component of index for key series
      amplee=0.5;
      lamdaspline = splinep(plfi,amplee); % spline parameter
      hkeeper=csaps(yrthish,rwin2,lamdaspline,yrthish); 
      H(i1H:i2H,scid)=hkeeper;
      clear hkeeper i1H i2H yrthish amplee lamdaspline;
      
     
      
   end;			% End of ksw2 while loop for first detrending
   
     
   % Update curve-fit flag
   fv(scid,:)='-f'; % change flag to indicate fit
   S(scid,1)=scid;  % put core sequence number in col 1 of S
   
   
   
   if strcmp(autocrank,'No');
      ksw1=menu('Select','Another Core ?','QUIT');
   else;
      disp(['Finished curve fitting core ' int2str(iauto) ' of ' int2str(ns)]);
      if iauto==ns; % have done all cores
         ksw1=2;
      end;
   end;
   close all;
end		% End for ksw1 while loop for working on this core

% Set flag for diff vs ratio
switch detmode;
case 1; % ratio
    S(:,6)=1;
case 2; % difference
    S(:,6)=2;
end;
    
% Update history 
c=clock;
c=num2str(c(4:5));
d=date;
Fwhen{3,2}=[d ', ' c];

% Save the vectors in a .mat file
nsv=menu('Save variables?','Add to original .mat storage file?',...
   'Store old data plus  new or revised G and S in a new file?','QUIT');
if nsv==1; % add variables to the original input .mat file
    Fwhen{3,4}=Fwhen{3,3};
   eval(['save ' pf1 ' Fwhen S  G ' ' -append']); % Fwhen is possibly new or revised
elseif nsv==2,
   [ofmat,path2]=uiputfile('*.mat','new .MAT file to store old data plus G, S, Fwhen: ');
   Fwhen{3,4}=ofmat; 
   pf2=[path2 ofmat];
   eval(['save ' pf2 ' X yrs nms A tstatus Fwhen S G ']);
end
% End of main function


% SUBFUNCTIONS

% subfn01 for individual series robust LFI
function [H,yrH,rath]=subfn01(X,yrs,S,G,A,coremask,plfi);
% H == time series matrix of low freq component of individual series (after detrending by current growth curve G)
% yrH == year vector for H
% X,yrs,G,A as in crvfit.m
% coremask 1 if use series, 0 if not.  If not, NaNs returned in col of H
% plfi period (yr) at which amp freq response of spline used for LFI is 0.5

nser=size(yrs,1); % number of series

% Compute size of tsm H and initialize
yrgo = min(yrs(:,1));
yrsp = max(yrs(:,2));
yrH = (yrgo:yrsp)';
mH = length(yrH);
H = repmat(NaN,mH,nser);
rath = repmat(NaN,nser,1); % to hold scaling ratio for robust LFI 

% Loop over series

for n = 1:nser;
    i1 = yrs(n,3);
    yron = yrs(n,1);
    yroff=yrs(n,2);
    yrthis = (yron:yroff)';
    nthis=length(yrthis);
    i2 = i1 + nthis-1;
    x = X(i1:i2); % the measured ring width
    g =G(i1:i2); % the previously fitted growth curve
    
    if ~isempty(A); % if want diff index
        xv=x;
        % Will be using matched transformed ring width instead of original ring width
        cshift = A{n,1}; % inc to add to rw before power tran
        ppower= A{n,2}; % power of transformation
        acoef = A{n,3}; % shift coef for matching
        bcoef = A{n,4}; % mult coef for matching
        
        % Power tran
        if ppower==1; % moot power; no transformation
            xvtran=xv+cshift;
        elseif ppower==0; % log transform
            xvtran=log10(xv+cshift);
        elseif ppower>0; % positive power
            xvtran= (xv+cshift) .^ ppower;
        elseif ppower<0; % neg power
            xvtran= -(xv+cshift) .^ ppower;
        end;
        % Matching
        xvmatch=acoef + bcoef * xvtran;
        xv=xvmatch;
        x=xv;
        y = x-g;
        depdiff = y - mean(y);  % as departures from mean
        depdiff = depdiff * S(n,7); % scaled departures
        y = depdiff +1.0; % with mean of 1.0
    else; % want ratio index
        y = x ./ g;
        % Adjust to exactly mean 1.0
        mny = nanmean(y);
        y = y-mny + 1.0;
    end;
    
    amp=0.5;
    p = splinep(plfi,amp); % spline parameter
    yorig=y;
    y=csaps(yrthis,y,p,yrthis); 
    rath(n) = std(y)/std(yorig); % ratio of std deviation of smoothed index to std  of original index
    
    
    % Put series into slot in tsm H
    L = yrH>=yron & yrH<=yroff;
    H(L,n)=y;
    
    
end;

disp(nser);





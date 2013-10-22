function crvfit1a
% crvfit1a: automatic fitting of growth curves to retain common low frequency variance
% crvfit1a
% Last revised 12-18-00
%
% Fits curves such that low frequency variance in longest series can be retained and used to guide fits
% for shorter series. 
%
%*** IN ****************************************
%
% User prompted for name of .mat file storing the ringwidth data
%    and associated years and names.  Ringwidth data previously put in
%       this file with rwlinp.m
% User prompted for several parameters
%   Number of straight line curves between horizontal line and least squares fit (default 100)
%   Most flexible spline to be entertained -- as fraction of series length (default 0.7 N)
%   Number of spline fits between the most flexible spline and the straight-line fit (defailt 100)
% User prompted to point to first master series (default is longest) 
% User prompted to click on optional .mat infile with order-of-treatment for curve fitting 
%   ( default longest to shortest)
%
%*** OUT ***********************************
%
% User prompted for name of .mat file storing curve-fit information. Can
% store in a new file, or add to the .mat file that holds the input
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
% cols 6-9 reserved for future use.  Previously these slots used to store
%  "second detrending" data analagous to cols 2-5.
% col 10, 11 -- start, end year of period for fitting curve
% col 12 -- starting row index in sov X for year in col 10
%
% Curves us full length of ring widths for fits
% 
% Sequence of steps
%
%   Store ring widths in a tsm
%   Compute and store overlap information
%   Compute candidate indices for all series, storing curve-fit params
%      C1 horizontal line  ( store mean level of ring width)
%      C2 straight lines (store slope and intercept)
%      C3 splines (store spline 'length' in years, fraction of N, spline parameter
%      The above are 3-d arrays, with i=series, j=param, and k=curve no. (C1 has just k=1)
%      D  indices; 3-d arrays with horiz line in bottom layer, straight lines above, splines above those
%      In D, i=series and j=year
%   Loop over master series in specified order, or from longest to shortest
%      Get all overlap series
%      Loop over candidate fits for master
%         Compute correlations between candidate master and all candidate ref curve fits
%         Find the candidate ref curve with highest r and store that rmax for each overlap series
%         Compute robust mean of rmax
%      Choose master curve-fit giving highest robust mean correlation
%   Store curve-fit info in form compatible with crvfit1b.m and other \STZ\ functions        



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


clear flmat



%***********************************************************************


%-- Prompt for all-the-same curve fits
autocrank=questdlg('Same curve-fit for all cores?');
if strcmp(autocrank,'Yes');
   kmen1=menu('Choose automatic curve fit for all series',...
      'Nonincreasing spline',...
      '0.8N spline',...
      'Horizontal Line');
   if kmen1==1;
      curvetype='sni';
   elseif kmen1==2; 
      curvetype='s8';
   elseif kmen1==3;
      curvetype='HL';
   end;
else;
   curvetype='null';
   kmen1=NaN;
end;

%-- Prompt for whether to use "blocking out".
if strcmp(autocrank,'No');
   blockuse=questdlg('Will you be blocking out segments in curve fitting?');
   blockuse=upper(blockuse); % change blockuse to upper case
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
   
   
   % Get ring-width series and vector of years; plot ring-width time series
   xv=X(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));
   gv=G(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1)); % growth trend, which
   %   might still be all NaN if no curve yet fit
   
   yrv=(yrs(scid,1):yrs(scid,2))';  % col vector of years for xv and gv
   
   if strcmp(autocrank,'No');
      hf0=figure('Units','normal','Position',[kfx,kfy,0.75,0.75]);
      if strcmp(fv(scid,:),'-n');
         plot(yrv,xv,'b');title(['Ring Width ',nms(scid,:)]);
      else;
         plot(yrv,xv,'b',yrv,gv,'m');title(['Ring Width ',nms(scid,:) '& curve fit']);
      end;
      
      xlabel('Year');ylabel('Ring Width, (0.01 mm)');
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
      else;
         ksw2=1;
         gv=repmat(NaN,length(gv),1);
      end;
   else; % autocrank mode
      ksw2=1;
      gv=repmat(NaN,length(gv),1);
   end;
   
      
         
   %**************  FIRST DETRENDING **************
   
   while ksw2==1;  % working on first detrending
      
      % Prompt for interval of series to fit detrending curve to
      if strcmp(autocrank,'Yes');
         nseg=4;
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
         if strcmp(curvetype,'s8') | strcmp(curvetype,'sni');
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
      
      % Store curve-type for detrending
      S(scid,2)=nfit1;
      
      if nfit1==9; % Cubic smoothing spline
         if strcmp(autocrank,'No');
            splp=menu('p-option ?','%N and .5 Amp','#yrs and .5 Amp','%N and x Amp',...
               '#yrs and x Amp','Specify p','Non-Increasing');
         else; % autocrank mode
            if strcmp(curvetype,'s8');
               splp=1;
            elseif strcmp(curvetype,'sni');
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
         else
            if scid==59;
               disp('here');
            end;
            
            if strcmp(autocrank,'No');
               [p,per]= monotspl(yrv,yrvn,xvn,length(xv),1);
            else;
               [p,per]= monotspl(yrv,yrvn,xvn,length(xv),2);

            end;
            
            amp = 0.5;
         end
         param(2:4)=[p per amp];
         cvx = (cfspl(p,yrv,yrvn,xvn))'; % Compute spline
         tstr1=['  SPL, p = ',num2str(p)];
         
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
      
      % Plot trend-line superposed on ring width
      if strcmp(autocrank,'No');
         figure(hf1);hold off;
         plot(yrvn,xvn,'b',yrvn,cvx(xind1:xind2),'k');
         % title(['Ring Width ',nms(scid,:),tstr1]);
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
            
            % Plot trend line superposed on ring width
            figure(hf1);hold off;
            plot(yrvn,xvn,'b',yrv(xind1:xind2),cvx(xind1:xind2),'k');
            title(['Ring Width ',nms(scid,:),tstr1]);
            xlabel('Year');ylabel('Ring Width, (0.01 mm)');
         end
         sgc=sgc+1;
      end 			% End of ksw3 while loop
      
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
            title([nms(scid,:),'- Tree-Ring Index by Ratio Method']);
            xlabel('Year');ylabel('Index');
         end;
         
      else; % want to re-do curve fit
         ksw2=1;
      end
      
      % Store the growth trend
      gv(xind1:xind2)=cvx(xind1:xind2);
      G(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1))=gv; %
      
      
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


% Save the vectors in a .mat file
nsv=menu('Save variables?','Add to original .mat storage file?',...
   'Store X, yrs, nms, G and S in a new file?','QUIT');
if nsv==1; % add variables to the original input .mat file
   eval(['save ' pf1 ' S  G ' ' -append']);
elseif nsv==2,
   [ofmat,path2]=uiputfile('*.mat','new .MAT file to store S,X,nms,yrs,G: ');
   pf2=[path2 ofmat];
   eval(['save ' pf2 ' X yrs nms S G ']);
end
% End of file

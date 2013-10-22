function crvfit
% crvfit: interactive fitting of curves to detrend ring-width series
% crvfit
% Last revised 9-2-99
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
% User defines endpoint of segment to be fit
% Detrending curve is fit to the segment
% User can block out sub-segments and re-fit smoothed curve to
% 	remaining data
% User accepts or revises detrending
%

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
if ~(exist('X')==1) | ~(exist('nms')==1) | ~(exist('yrs')==1),
  error('The selected .mat file does not contain X, nms and yrs');
end
clear flmat


%-- Prompt for whether to use "blocking out".
blockuse=questdlg('Will you be blocking out segments in curve fitting?');
blockuse=upper(blockuse); % change blockuse to upper case

%***********************************************************************


% Get the number of cores, which equals the row size of  matrix
% of core ids
ns=size(nms,1);

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

while ksw1==1;  % while working on this core
   nsms=[nsmt nms fv]; % string matrix with sequence number,
   % nsms give core id, and fit-status: -f == "fit", '-n'="not yet fit" 
   
   % Select a core.  
   scid=menu('Core ID # ?',(cellstr(nsms))');  % sequence number of selected core
   
   % Get ring-width series and vector of years; plot ring-width time series
   xv=X(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));
   yrv=(yrs(scid,1):yrs(scid,2))';  % col vector of years for xv
   hf0=figure('Units','normal','Position',[kfx,kfy,0.75,0.75]);
   plot(yrv,xv,'b');title(['Ring Width ',nms(scid,:)]);
   xlabel('Year');ylabel('Ring Width, (0.01 mm)');
   
   xv1=xv; % store full-length ring-width series
   
   
   ksw2=1;  % while flag for trial and  detrending 
      
   %**************  FIRST DETRENDING **************
   
   while ksw2==1;  % working on first detrending
      
      % Prompt for interval of series to fit detrending curve to
      nseg=menu('Ends ?','Graphical input','Specify in #','Full Length');
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
      hf1=figure('Units','normal','Position',[kfx,kfy,0.75,0.75]);
      plot(yrvn,xvn,'b');title(['RingWidth ',nms(scid,:)]);
      xlabel('Year');ylabel('RW, (0.01 mm)');
      
      param=repmat(NaN,1,8);  % Initialize parameters for curve fit
      
      % Prompt for curve fit options and Curve fitting
      % eventual curve options ['1';'2';'3';'4';'5';'6';'7';'8';'9'];
      nopt1=['1 - Neg Ex';'2 - SL    ';'4 - Mean  ';'9 - Spline'];
      
      % Set flag for type of detrending curve
      nfit = menu('Curve-fit option',(cellstr(nopt1))');

      %nfit=svmenu('1st fit option',nopt1);
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
         splp=menu('p-option ?','%N and .5 Amp','#yrs and .5 Amp','%N and x Amp',...
            '#yrs and x Amp','Specify p','Non-Increasing');
         if splp>7;
            jh1=jdisp('Please return to command window for input');
            pause(1);
            close(jh1);
         end;
                
         
         if splp==1; % you specify the wavelength with 0.5 AFR as fraction of series length
            prompt={'Enter the decimal fraction of series length:'};
            def={'.70'};
            dlgTitle='0.50 Frequency Response Wavelength';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            pper = str2num(answer{1});
            %pper=input('Period decimal fraction of the total series length = ');
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
            [p,per]= monotspl(yrv,yrvn,xvn,length(xv));
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
      
      % Plot first trend-line superposed on ring width
      figure(hf1);hold off;
      plot(yrvn,xvn,'b',yrvn,cvx(xind1:xind2),'k');
      % title(['Ring Width ',nms(scid,:),tstr1]);
      xlabel('Year');ylabel('Ring Width, (0.01 mm)');
      
      % Initialize blocking settings
      ksw3=1; % while control for blocking segments from use in curve fit
      nblk=1;
      sgc=1;
      
      %********** WHILE BLOCKING OUT INTERVALS ***********
      
      while ksw3==1 & nblk~=6 & ~strcmp(blockuse,'YES');
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
      ksw2=menu('Select One','Change curve-fit ?','Accept curve-fit'); 
      if ksw2==2,
         nq=logical(1);
      elseif ksw2==1;
         nq=logical(0);
      else
         error('Impossible value for ksw2');
      end
      
      if nq; % you accepted curve-fit
         % compute index as ratio of ring-width to fitted curve
         rwin1=xv1(xind1:xind2)./cvx(xind1:xind2);
         % Comnpute  mean of ratio ring-width to fitted curve 
         mnx=nanmean(rwin1);
         % compute index adjusted to mean 1.0
         rwin2 = rwin1 - mnx +1.0;
         
         % Plot index
         hf2=figure('Units','normal','Position',[kfx,kfy,0.75,0.75]);
         plot(yrv(xind1:xind2),rwin2,'b',[yrv(xind1) yrv(xind2)],[1.0 1.0],'r');
         title([nms(scid,:),'- Tree-Ring Index by Ratio Method']);
         xlabel('Year');ylabel('Index');
      else; % want to re-do curve fit
         ksw2=1;
      end
   end;			% End of ksw2 while loop for first detrending
   
     
   % Update curve-fit flag
   fv(scid,:)='-f'; % change flag to indicate fit
   S(scid,1)=scid;  % put core sequence number in col 1 of S
   
   ksw1=menu('Select','Another Core ?','QUIT');
   close all;
end		% End for ksw1 while loop for working on this core


% Save the vectors in a .mat file
nsv=menu('Save variables?','Add to original .mat storage file?',...
   'Store X, yrs, nms, and S in a new file?','QUIT');
if nsv==1; % add variables to the original input .mat file
   eval(['save ' pf1 ' S ' '-append']);
elseif nsv==2,
   [ofmat,path2]=uiputfile('*.mat','new .MAT file to store S,X,nms,yrs: ');
   pf2=[path2 ofmat];
   eval(['save ' pf2 ' X yrs nms S']);
end
% End of file

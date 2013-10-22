function cfit1a

% DESCRIPTION : cfit1
% Interactive fitting of detrending curves by double detrending
% Curve-fit info will be stored in S
%
% INPUT  :  NONE
% OUTPUT :  NONE
%
%
% Best clear the workspace before starting.
% User defines endpoint of segment to be fit
% First-detrending curve is fit to the segment
% User can block out sub-segments and re-fit smoothed curve to
% 	remaining data
% User accepts first detrending
% Choose curve-type for second detrending
% Second detrend, assuming same blocked sub-segments
% 
%_____________________________________________________________

notan = NaN; 

clear all ; % clear all variables from the workspace

% Prompt for name of .mat file with ring widths, core ids, and 
% year information
flmat=uigetfile('*.mat','.MAT filename ?');

% Do strung-out ring-width vector X, and yrs and nms exist in the workspace?
if ~exist('X') | ~exist('nms') | ~exist('yrs'),
  eval(['load ',flmat]);
end

% Get the number of cores, which equals the number of rows in the matrix
% of core ids
[ns,dum1]=size(nms);  %  ns is number of cores


% Get information on whether particular cores have already been fit
% in a previous run.  If this is a first run, initialize S and the
% string matrix fv. If S exists, col 1 of S will either be a sequence
% number or 0, depending on whether the core has been fit or not.

fvv='-n'; %  not fit yet
fvf='-f'; %  already fit 

% If S exists, get the info on whether cores have been fit and store in fv,
% If S does not exist, initialize as zeros
fv = fvv(ones(ns,1),:); % initialize string matrix as if no cores yet fit
if exist('S'),
  LS1=S(:,1)~=0;  % non-zero in first col of S:  core had been fit
  numf = sum(LS1);  % number of cores already fit
  fv (LS1,:) = fvf(ones(numf,1),:);
else; % Intialize matrix S as zeros
  S=zeros(ns,22);
end



% Make a string matrix with the following for each row:
%  * the sequential number of the core
%  * '-'
%  So, might get
%	1-
% 	2-
%  etc
nsmt=[];
for i=1:ns,
 nsmt=str2mat(nsmt,[int2str(i),'-']);
end
nsmt(1,:)=[];  % To remove row 1, which is blanks



ksw1=1; % while switch for going to another core

kfx=0.2; % x and y plotting positions for tsp of ringwidth and fitted curve
kfy=0.15;


%******************  WHILE OVER CORES  ********************

while ksw1==1;  % while working on this core

  nsms=[nsmt nms fv]; % string matrix with sequence number,
	% core ids, and fit-status of cores
 % Seclect a core
 scid=svmenu('Core ID # ?',nsms);
 
 
 % Get ring-width series and vector of years; plot ring width time series
 xv=X(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));
 yrv=(yrs(scid,1):yrs(scid,2))';  % col vector of years for xv
 hf0=figure('Color','w','Units','normal','Position',[kfx,kfy,0.75,0.75]);
 plot(yrv,xv,'b');title(['Ring Width ',nms(scid,:)]);
 xlabel('Year');ylabel('Ring Width, (0.01 mm)');


 xv1=xv; % full-length ring-width series

 ksw2=1;  % Switch for first detrending in double-detrending seequence


%************** WHILE FIRST DETRENDING **************

 while ksw2==1;  

  % Prompt for interval of series to fit
  nseg=shlmnu('Ends ?','Graphical input','Specify in #','Full Length');
  if nseg==1,
    jh1=jdisp('Please click at two corner points');
    pause(1);
    close(jh1);
    figure(hf0);
    [segx,segy]=ginput(2);
    [erflg,xind1,xind2]=erchk(yrv,min(segx),max(segx));
  elseif nseg==2,
    jh1=jdisp('Please return to command window for input');
    pause(1);
    close(jh1);
    yeargo=input('Please enter begin year = ');
    yearstop=input('Please enter end year = ');
    [erflg,xind1,xind2]=erchk(yrv,yeargo,yearstop);
  else
    xind1=1;
    xind2=length(yrv);
    erflg=0;
  end

  if erflg==-1,
    jh1=jdisp('Years out of range. Re-run the program');
    pause(1);
    close(jh1);
    return;
  end

  yrvn=yrv(xind1:xind2); % year vector for selected interval
  xvn=xv(xind1:xind2);   % ring-width data for the selected interval

  S(scid,10)=yrv(xind1); % Store start year for interval
  S(scid,11)=yrv(xind2); % Store end year for interval
  S(scid,12)=yrs(scid,3)+xind1-1; % store starting index for interval

  % Plot time series of ring width for selected segment
  hf1=figure('Color','w','Units','normal','Position',[kfx,kfy,0.75,0.75]);
  plot(yrvn,xvn,'b');title(['RingWidth ',nms(scid,:)]);
  xlabel('Year');ylabel('RW, (0.01 mm)');

  param=notan(:,ones(8,1));  % Initialize parameter matrix whose values
		% will be stored in S

  % Prompt for curve fit options and Curve fitting
  % eventual curve options ['1';'2';'3';'4';'5';'6';'7';'8';'9'];
  nopt1=['1 - Neg Ex';'4 - Mean  ';'9 - Spline'];
  
  nfit=svmenu('1st fit option',nopt1);
  if nfit==1; % neg exp
    nfit1=1;
  elseif nfit==2; % horiz thru mean
    nfit1=4;
  else
    nfit1=9; % spline
  end
  S(scid,2)=nfit1;

  if nfit1==9; % Cubic smoothing spline
    splp=shlmnu('p-option ?','%N and .5 Amp','#yrs and .5 Amp','%N and x Amp',...
                  '#yrs and x Amp','Specify p','Non-Increasing');
    jh1=jdisp('Please return to command window for input');
    pause(1);
    close(jh1);
    if splp==1,
      pper=input('Period as % of the total series length = ');
      per=pper*length(xv);
      amp=0.5;
      p=splinep(per,amp);
    elseif splp==2,
      per=input('Length of the period (Years) = ');
      amp=0.5;
      p=splinep(per,amp);
    elseif splp==3,
      pper=input('Period as % of the total series length = ');
      per=pper*length(xv);
      amp=input('Please enter the value of Amplitude = ');
      p=splinep(per,amp);
    elseif splp==4,
      per=input('Length of the period (Years) = ');
      amp=input('Please enter the value of Amplitude = ');
      p=splinep(per,amp);
    elseif splp == 5,
      p=input('Spline parameter, p = ');
      per=NaN;
      amp=NaN;
    else
      [p,per]= monotspl(yrv,yrvn,xvn,length(xv));
      amp = 0.5;
    end
    param(2:4)=[p per amp];
    cvx = cfspl(p,yrv,yrvn,xvn); % Compute spline
    tstr1=['  SPL, p = ',num2str(p)];

  elseif nfit1==1; % Neg exp
    
    hwarn = 1;  % want to display warning dialog if neg exp wrong type
    [cvx,k,a,b] = cfnegx(yrv,yrvn,xvn,hwarn);
    param(2:4)=[k a b];
    tstr1 = '  NEG EXP ';

  elseif nfit1==4; % Horizontal line through mean
    cvx = cfmean(yrv,yrvn,xvn); % Horizontal line through mean
    tstr1='  HMN ';
    param(2)=mean(xvn(~isnan(xvn))); % sample mean of the valid ringwidths
  end

  S(scid,3:5) = param(2:4); % Store parameters
  
  % Plot first trend-line superposed on ring width
  figure(hf1);hold off;
  plot(yrvn,xvn,'b',yrvn,cvx(xind1:xind2),'r');
  title(['RingWidth ',nms(scid,:),tstr1]);
  xlabel('Year');ylabel('RW, (0.01 mm)');



  % Initialize blocking settings
  ksw3=1;  
  nblk=1;
  sgc=1;

%********** WHILE BLOCKING OUT INTERVALS ***********

  while ksw3==1 & nblk~=6,
   % Blocking out specified data segments
   if sgc==1,
     xv1=xv;
   end
   [S,yrvn,xvn,xv1,nblk] = blocdat(hf1,scid,xind1,yrvn,...
     xvn,yrv,xv,xv1,sgc,S);
   if nblk==6, break; end
   figure(hf1);hold off;
   plot(yrvn,xvn,'b');title(['RingWidth ',nms(scid,:)]);
   xlabel('Year');ylabel('RW, (0.01 mm)');
   pause(1);
   ksw3=shlmnu('Select','Block new segments ?','Continue');

   if ksw3~=1,
    % Curve fitting for the blocked data
    if nfit1==9; % Spline
       cvx = cfspl(p,yrv,yrvn,xvn);
    elseif nfit1==1; % Neg exp
      hwarn = 1;  % want to display warning dialog if neg exp wrong type
      [cvx,k,a,b] = cfnegx(yrv,yrvn,xvn,hwarn);
      param(2:4)=[k a b];
      S(scid,3:5) = param(2:4); % Store parameters
      tstr1 = '  NEG EXP ';
    elseif nfit1==4; % Horizontal line through mean
       cvx = cfmean(yrv,yrvn,xvn);
    end

    % Plot first-detrending line superposed on ring width
    figure(hf1);hold off;
    plot(yrvn,xvn,'b',yrv(xind1:xind2),cvx(xind1:xind2),'r');
    title(['Ring Width ',nms(scid,:),tstr1]);
    xlabel('Year');ylabel('Ring Width, (0.01 mm)');
   end
   sgc=sgc+1;
  end 			% End of ksw3 while loop


  % If 
  ksw2=shlmnu('Select','Change Options ?','Continue');
  if ksw2~=1,
   % Accept the current 1st detrending ?
   nq=usinp('Accept 1st detrending ?');
   if nq,
    rwin1=xv1(xind1:xind2)./cvx(xind1:xind2);
    hf2=figure('Color','w','Units','normal','Position',[kfx,kfy,0.75,0.75]);
    plot(yrv(xind1:xind2),rwin1,'r');
    title([nms(scid,:),'- Index After First Detrending']);
    xlabel('Year');ylabel('Index');
   else
    ksw2=1;
   end
  end
 end			% End of ksw2 while loop


%****************  SECOND DETRENDING ********************

% If first detrending was by negative exponential, this section 
% allows interactive second detrending by spline (ratio method), or
% subtraction of the mean. If first detrending was not by a negative
% exponential, this "second detrending" amounts merely to subtraction
% of the sample mean followed by addition of 1.0 -- the final 
% index therefore has a mean of 1.0 over the full length of the 
% core's ring-width series. The second-curve specification of "4" means
% this additive conversion to mean of 1, and is done automatically
% without interaction for first-curve fits other than the negative 
% exp.

 
 nopt2=['4 - Mean  ';'9 - Spline'];


 if nfit1~=1; %  if not neg exp first fit, convert to zero mean additively
  xvn=rwin1;
  param(1:4)=[notan(:,ones(4,1))]; % Initialize param(1:4) as NaNs
  LL = isnan(xvn);
  LL = ~LL;
  xgood = xvn(LL);
  nfit2=4; % subtract mean
  rwin2=rwin1-mean(xgood)+1.0;
  fv(scid,:)='-f'; % change flag to indicate fit
  S(scid,1)=scid;  % put core sequence number in col 1 of S
  S(scid,6)=4;
  S(scid,7)=mean(xgood);
 else;  % neg exp was used for first fit

   nq=0;
   while nq~=1,
    nfit=svmenu('2nd fit option',nopt2);
    if nfit==1,
     nfit2=4;
    elseif nfit==2,
     nfit2=9;
    end

    param(1:4)=notan(:,ones(4,1));
    xvn=rwin1; % Will detrend the first index

    if nfit2==9;  % Spline
      splp=shlmnu('p-option ?','%N and .5 Amp','#yrs and .5 Amp','%N and x Amp',...
                  '#yrs and x Amp','Specify p');
      jh1=jdisp('Please return to command window for input');
      pause(1);
      close(jh1);
      if splp==1,
        pper=input('Period as % of the total series length = ');
        per=pper*length(xv);
        amp=0.5;
        p=splinep(per,amp);
      elseif splp==2,
        per=input('Length of the period (Years) = ');
        amp=0.5;
        p=splinep(per,amp);
      elseif splp==3,
        pper=input('Period as % of the total series length = ');
        per=pper*length(xv);
        amp=input('Please enter the value of Amplitude = ');
        p=splinep(per,amp);
      elseif splp==4,
        per=input('Length of the period (Years) = ');
        amp=input('Please enter the value of Amplitude = ');
        p=splinep(per,amp);
      else
        p=input('Spline parameter, p = ');
        per=NaN;
        amp=NaN;
      end
      param(2:4)= [p per amp];
      cvx2 = cfspl(p,yrv,yrvn,xvn);
      tstr2=['  SPL, p = ',num2str(p)];

    elseif nfit2==4; % Horizontal line through mean
      xgood =xvn(~isnan(xvn));
      rwin2=rwin1-mean(xgood)+1.0;
      param(2)=mean(xgood);
    end

    S(scid,6)=nfit2;
    S(scid,7:9)=param(2:4);
  
    if nfit2==9,
      hf3=figure('Color','w','Units','normal',...
	  'Position',[kfx,kfy,0.75,0.75]);
      plot(yrvn,xvn,'r',yrv(xind1:xind2),cvx2(xind1:xind2),'b');
      title(['2nd curve fit ',nms(scid,:),tstr2]);
      xlabel('Year');
      % Accept the current 2nd detrending ?
      nq=usinp('Accept 2nd detrending ?');
    elseif nfit2==4
      nq=1;
    end


    if nq,
     fv(scid,:)='-f'; % change flag to indicate fit
     S(scid,1)=scid;  % put core sequence number in col 1 of S
     
    
     if nfit2==9,
       hf4=figure('Color','w','Units','normal','Position',[kfx,kfy,0.75,0.75]);
       rwin2=rwin1./cvx2(xind1:xind2);
       plot(yrv(xind1:xind2),rwin2,'r'),
       title([nms(scid,:),'- Index after 2nd fit']);
       xlabel('Year');ylabel('Index');
     end
    end
   end            % End of nq while loop; accepted second detrending
 end  % of if nfit~=1
 ksw1=shlmnu('Select','Another Core ?','QUIT');
 close all;
end		% End for ksw1 while loop

nsvf=usinp('Save data in file ?');
if nsvf,
  flsn=uiputfile('*.mat','Save as ?');
  eval(['save ' flsn ' X ' ' yrs ' ' nms ' ' S ']);
end


% End of file

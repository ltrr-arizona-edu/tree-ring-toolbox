function pdgm5(x,yr,uiinf)
% pdgm5:  smoothed periodogram spectral analysis, univariate
% CALL: pdgm5(x,yr,uiinf);
%
% Meko 8-26-98
%
%******************  INPUT
%
% x (mx x 1)r   time series
% yr (mx x 1)r   year vector for x
%
% uiinf{}  user-interface input information that can either be passed as
%   input argument, or can be prompted for by pdgm5.m.   If uiinf is [],
%   pdgm5 will prompt you for the information.  Otherwise, the contents
%   of uiinf{} should be as follows:
%
%  {1} name (1 x ?)s   name of series
%  {2} yrsa (1 x 2)i  start and end years of analysis period
%  {3} pdvar (1 x 2)r  longest and shortest period lengths (yr) bracketing
%   period range want percentage of variance computed for (e.g., [inf 5])
%  {4} dan1 (1 x ?)i widths of daniel filters to be used for spectral estimates
%  {5} dan2 (1 x ?)i  widths of daniel filters to be used for null continuum
%  {6} kopt (1 x 1)i  options
%    kopt(1): hypothesis test mode"
%       ==1 yes
%       ==0 no


%############### DESCRIPTION ##########################################
%
% Method follows procedures described in the book "Fourier Analysis of 
% Time Series:  An Introduction", by Peter Bloomfield, 1976, John Wiley
% & Sons
%
% The mains steps are as follows:
%
% 1. Read in, subtract the mean from, taper, and pad a time series.
% 2. Compute the discrete Fourier transform (DFT) of the series.
% 3. Compute the periodogram from the DFT.
% 4. Smooth the periodogram with a series of Daniell filters to get
%	the estimated spectrum.  
% 5. Smooth the periodogram with a separate set of (broader) Daniell
%	filters to get a smooth null continuum.
% 6. Compute 95% confidence interval around the estimated spectrum using
%	the chi-squared distribution.
% 7. Summarize the spectral peaks in table form.
%


%*************** GRAPHICS WINDOWS
%
windlist = 'Windows';
windlist=char(windlist,...
   '1 - tsp of raw time series, analysis period',...
   '2 - tsp of raw time series with mean subtracted',...
   '3 - data taper window',...
   '4 - tsp of tapered, padded data',...
   '5 - discrete Fourier Transform',...
   '6 - raw (before smoothing) periodogram',...
   '7 - current version of smoothed periodogram (spectrum) overlaid on raw pdgm',...
   '8 - current version of null continuum overlaind on raw pdgm',...
   '9 - final version of spectrum, with 95% CI and null continuum',...
   '10 - text summary window');


%************************* Size time series
[m1,n1]=size(x);
[m2,n2]=size(yr);
if n1~=1 | n2~=1 | m1~=m2 | any(isnan(x));
   error('x and yr must be col vectors of same length, and no NaN in x');
end
nyr = length(x);
yrs = [min(yr)  max(yr)];



%************************  GET INPUT CONTROL (ALL EXCEPT THE DANIEL FILTER LENGTHS)

%--------  Unload uiinf, if uiinf is not empty
if ~isempty(uiinf);
   kmode=2;
   name = uiinf{1}; % name of time series (used in figure windows)
   yrsa = uiinf{2}; % period for analysis
   nyra=yrsa(2)-yrsa(1)+1;
   pdvar = uiinf{3}; % perio range for variance computation
   dan1 = uiinf{4};  % daniel widths for smoothing pdgm into spectrum
   dan2 = uiinf{5}; % daniel widths for smoothing pdgm into null continuum
   kopt = uiinf{6}; % options
else
   kmode=1;
   % series name, analysis period, period for variance computation
   prompt={'Enter the series name:','Enter the period for analysis',...
         'Enter the period range (yr) for variance computation:'};
   def={'Series 1',int2str(yrs),'[inf 4]'};
   titpmt1='Input Information';
   lineNo=1;
   answer=inputdlg(prompt,titpmt1,lineNo,def);
   name = answer{1};
   yrsa = str2num(answer{2});
   nyra=yrsa(2)-yrsa(1)+1;
   pdvar = str2num(answer{3});
      
   % Daniel filter lengths
   prompt={'For smoothing pdgm into spectral estimates:',...
         'For smoothing pdgm into null continuum'};
   % compute default for null continuum daniel filter lengths, ensuring that they odd
   d1 = round(nyra/7); d2=round(nyra/5);
   if mod(d1,2)==0;
      d1=d1-1;
   end
   if mod(d2,2)==0;
      d2=d2-1;
   end
   def={'[7 11]',['[' int2str([d1 d2]) ']']};
   titpmt2='Input the Daniel Filter Lengths';
   lineNo=1;
   answer=inputdlg(prompt,titpmt2,lineNo,def);
   dan1= str2num(answer{1});
   dan2 = str2num(answer{2});
    
end

if ~ischar(name);
   error('Series name not character');
end
if yrsa(1)<min(yr) | yrsa(2)>max(yr);
   error('Specified period for analyis outside time range of time series');
end
yra = (yrsa(1):yrsa(2))'; % year vector for analysis period
if pdvar(1)<=pdvar(2);
   error('pdvar should have longer period as first element');
end
if any(pdvar<2);
   error('pdvar has a value less than 2 yr');
end

if any(mod(dan1,2)==0) | any(mod(dan2,2)==0);
   error('Lengths of Daniel filters must be odd');
end

% Hypothesis test mode:
khypo = questdlg('Hypothesis Test Mode?');
if strcmp(khypo,'Yes');
   kmen9 = menu('Frame Hypothesis in Terms of What?',...
      'Period','Frequency');
   if kmen9==1;
      prompt={'H0: Spectrum no different than null continuum at this period (yr):'};
      def={'20'};
      tit='Null Hypothesis';
      lineNo=1;
      answer=inputdlg(prompt,tit,lineNo,def);
      keypd=str2num(answer{1});
      keyfreq = 1/keypd; % frequency (cycles per year)
   elseif kmen9==2;
      prompt={'H0: Spectrum no different than null continuum at this frequency (1/yr):'};
      def={'.05'};
      tit='Null Hypothesis';
      lineNo=1;
      answer=inputdlg(prompt,tit,lineNo,def);
      keyfreq=str2num(answer{1});
      keypd = 1/keyfreq; % frequency (cycles per year)
   end
else
   keypd=[];
   keyfreq=[];
end

clear d1 d2 answer titpmt1 titpmt2 m1 m2 n1 n2 lineNo prompt def;



%******************* GET SEGMENT OF TIME SERIES FOR ANALYSIS

L1 = yr>=yrsa(1) & yr<=yrsa(2);  % pointer to analysis period in x,yr
x = x(L1);
yr = yr(L1);

stra1=sprintf('%5.0f-%5.0f',yrsa);
stra1 = ['Analysis Period: ' stra1];
stra2 = ['\itN = \rm' int2str(nyra) ' yr'];


%***************** PROMPT FOR HOW MUCH OF ENDS OF SERIES TO TAPER, AND FOR PADDED LENGTH

% taper
prompt={'Enter total decimal fraction of series to taper: '};
def={'.10'};
tit='Taper Fraction';
lineNo=1;
answer=inputdlg(prompt,tit,lineNo,def);
mtap1=str2num(answer{1});

% Padding
prompt={'Enter desired padded length (yr): '};
def={num2str(2^nextpow2(diff(yrsa)+1))};
tit='Desired Padded Length';
lineNo=1;
answer=inputdlg(prompt,tit,lineNo,def);
padlen=str2num(answer{1});


%*******************  PLOT OF RAW TIME SERIES

% Compute mean, standard dev, and variance of time series
xmn = mean(x);
xstd = std(x);
xvar = var(x); % variance of time series, before tapering

figure(1);
hp1=plot(yr,x);
cord=get(gca,'ColorOrder');
title(['Time Series Plot (' stra1 ')']);
xlabel('Year');
ylabel('Value'); 
xlim=get(gca,'XLim');
hl1=line(xlim,[xmn xmn]);
set(hl1,'color',cord(2,:));
legend(name,'Mean');
grid;
zoom xon;

%*******************  PLOT OF DETRENDED (MEAN SUBTRACTED) TIME SERIES

x=x-xmn; % detrended

figure(2);
hp1=plot(yr,x);
title([name ',  With Mean Subtracted (' stra1 ')']);
xlabel('Year');
ylabel('Value'); 
xlim=get(gca,'XLim');
hl1=line(xlim,[0  0]);
set(hl1,'color',cord(2,:));
grid;
zoom xon;


%****************** PLOT OF DATA TAPER WINDOW
figure(3);

% Taper mtap1/2 of data on each end.  See Bloomfield, p. 84.
% Compute data window
ends= fix((mtap1/2)*nyra);
dwin=ones(1,nyra);
t=1:ends;
btemp=0.5*(1-cos((pi*(t-0.5))/ends));  %intermediate variable
t=nyra-ends+1:nyra;
ctemp=0.5*(1-cos((pi*(nyra-t+0.5))/ends));  % intermediate

dwin(1:ends)=btemp;
dwin(nyra-ends+1:nyra)=ctemp;

hp3=plot(yr,dwin);
title(['Data Taper Window, ' stra1])
set(gca,'YLim',[0 1.1]);
set(hp3,'Linewidth',1.5);
xlabel('Year');
ylabel('Weight');
grid;
zoom xon;


%***************** TSP OF TAPERED DATA
figure(4);

% Taper by applying data window to data
x1=x .* dwin';
hp4=plot(yra,x1);
title(['Tapered Time Series ('  stra1 ')']);
xlabel('Year');
ylabel('Departure from Mean');
grid ;
zoom xon;

% Compute variance of tapered series
xvartap=var(x1); % 


%***************** PAD TIME SERIES WITH ZEROS ON RIGHT HAND SIDE
figure (5);

atemp=padlen-nyra;
btemp=zeros(atemp,1);
x2=[x1;btemp]; % padded, tapered time series for analysis period

yrsb = [yrsa(1) yrsa(2)+atemp];
yrb=(yrsb(1):yrsb(2))';

plot(yrb,x2);
title(['Tapered, Padded Time Series, ' name]);
xlabel('Year');
ylabel('Departure from Mean');
zoom xon;
grid;

stryr1 = sprintf('%4.0f',nyra);
stryr2 = sprintf('%4.0f',padlen);
stryr=['\itN=\rm ' stryr1 ',  \itN_{pad}=\rm ' stryr2];

%******************* PLOT DISCRETE FFT
figure(6);

% Compute FFT
% Note, what follows is same as in  Ljung Notation
% Bloomfield's squared dft is 1/N times Ljung's DFT.


z=fft(x2,padlen);

t=0:1:padlen;
f=t/padlen;

plot(f(1:padlen/2),real(z(1:padlen/2)));
title(['Real Part of Discrete Fourier Tansform of ' name]); 
xlabel('Frequency (yr^{-1})');



%***************** COMPUTE AND PLOT RAW PERIODOGRAM
figure(7);

% Compute periodogram by squaring DFT of tapered, padded series.
% Computation as defined by Ljung by squaring the DFT
% Note: Bloomfield's periodogram is 1/2*pi times Ljung's periodogram.

Pyy = z.*conj(z)/padlen;
hp7=plot(f(1:padlen/2),Pyy(1:padlen/2),'+g');
set(hp7(1),'Color',cord(3,:));
title(['Raw, or Unsmoothed, Periodogram, ' name]);
xlabel('Frequency (yr^{-1})');

% for convenience, store periodogram (first half) of values in y, and 
% associated frequencies in ff
ff=f(1:padlen/2);
y= Pyy(1:padlen/2);  % put pdgm ordinates ...


%*************** VARIANCE CHECK: COMPARE COMPUTED VARIANCE WITH PDGM SUM
clc;
xvarpdgm = (1/(nyra-1)) * sum(Pyy(1:padlen));
strv1 = sprintf('Variance computed directly from (tapered)time series = %10g',xvartap);
strv2 = sprintf('Variance computed from periodogram = %10g',xvarpdgm);
clear Pyy f


%**************** COMPUTE PERCENTAGE OF VARIANCE IN A PARTICULAR PERIOD RANGE
p1=pdvar(1); % longest period (yr)
p2=pdvar(2); % shortest period (yr)
f1=1.0/p1;    % frequency corresp to longest wavelength of range
f2=1.0/p2;    % frequency corresp to shortest ...
atemp=ones(padlen/2,1);

i1=find(ff>=f1&ff<=f2);  %  subscripts of ff in desired freq range
[m,n]=size(i1);          % n is number of elements found

pct=sum(y(i1))/sum(y);
%revise the "end periods" to agree with frequencies that have pdgm estimates at
if ff(i1(1))==0;
   p1=inf;
else
   p1=1.0/ff(i1(1));
end
if ff(i1(n))==0;
   p2=inf;
else
   p2=1.0/ff(i1(n));
end
strvpct1=sprintf('%5.1f pct of variance',pct*100);
strvpct2=sprintf('%7.2f yr -%7.2f yr',p1,p2);
strvpct=[strvpct1 ' in period range ' strvpct2];

%*************** EXTEND PERIODOGRAM BEYOND FREQUENCIES 0 AND 0.5 BEFORE SMOOTHING

%  Before smoothing periodogram, want to extend it so that it is 
% symmetric about frequencies 0 and 0.5  /yr
n=length(y);		% first element holds no ordinates in half-pdgm
yrev=y(n(1):-1:1);
yext=[yrev;y;yrev];
%plot(yext)
%title('Extended periodogram, symmetric about ends')
% A symmetric series of length n has been tacked on each end of the pdgm.
% So the original periodogram really covers the elements n+1 through 
% 2n of yext.  And, when convolute with an m-weight Daniell filter, will
% have effect of shifting (m-1)/2 units relative to central weight.  So will
% eventually want to plot and analyze this range of yext as the smoothed 
% periodogram:   n+1+(m-1)/2  thru  2n+(m-1)/2

% Always pick m to have odd number of weights.



%**********************  SMOOTHED PERIODOGRAM

ksfirst=1;
kwh1 = 1;
while kwh1==1;
   kmen1=menu('Choose 1','Smooth Periodogram or Revise Smoothing',...
      'Accept Spectrum');
   if kmen1==1;
      figure(7);
      if ksfirst==1;
         daniel1=dan1;
         ksfirst=2;
      else
         prompt={'For smoothing pdgm into spectral estimates:'};
         def={['[' int2str(daniel1) ']']};
         tit1='Input the Daniel Filter Lengths';
         lineNo=1;
         answer=inputdlg(prompt,tit1,lineNo,def);
         daniel1= str2num(answer{1});
      end
      w = danwgtn(daniel1);
      mm=length(w); % span of daniel filter
      % Compute bandwidth and df of the Daniell filter
      [bw,df]=danbw2(w,nyra,padlen);
      igo = n+1+(mm-1)/2;
      istop=2*n+(mm-1)/2;
      ysm=conv(yext,w); % smoothed periodogram == spectrum
      hp7b=plot(ff,y,'.g',ff,ysm(igo:istop),'-b');
      set(hp7b(1),'Color',cord(3,:)); % periodogram
      set(hp7b(2),'Color',cord(1,:)); %  spectrum
      
      % Plot bandwidth bar.  Change next 2 lines as needed for formatting
      xlim7=get(gca,'XLim'); 
      ylim7=get(gca,'YLim');
      delx = diff(xlim7);
      dely = diff(ylim7);
      ybwgo = ylim7(1)+0.83*dely;  % bar 8/10 of way to top of figure
      yinc = dely/150; % vertical width of tick at ends ob bw bar
      xbwgo=xlim7(1)+delx/2 - bw/2;
      xbwsp = xbwgo+bw;
      hlinebar1=line([xbwgo xbwsp],[ybwgo ybwgo]); % the bw bar itself (horiz line)
      hlinebar2=line([xbwgo xbwgo],[ybwgo-yinc/2 ybwgo+yinc/2]); % left side vert line
      hlinebar3=line([xbwsp xbwsp],[ybwgo-yinc/2 ybwgo+yinc/2]); % right side vert line
      set(hlinebar1,'Color',1-get(gca,'Color'));
      set(hlinebar2,'Color',1-get(gca,'Color'));
      set(hlinebar3,'Color',1-get(gca,'Color'));
      
      text(xbwgo+(xbwsp-xbwgo)/2,ybwgo+yinc/2,'BW',...
         'VerticalAlignment','bottom','HorizontalAlignment','center');
      
      xlabel('Frequency (yr^{-1})');
      ylabel('Relative Variance');
      grid;
      set(gca,'Position',[0.1300    0.1100    0.7750    0.78]);
      title({['Smoothed Pdgm of ' name ';  ' stryr],...
            ['Spans: [' int2str(daniel1) ']']});
      set(hp7b(2),'LineWidth',1.5);
   elseif kmen1==2;
      kwh1=0;
      ysm=ysm(igo:istop);
   end
end


%**********************  NULL CONTINUUM

ksfirst=1;
kwh1 = 1;
while kwh1==1;
   kmen1=menu('Choose 1','Make or Revise Null Continuum',...
      'Accept Null Continuum');
   if kmen1==1;
      figure(7);
      if ksfirst==1;
         daniel2=dan2;
         ksfirst=2;
      else
         prompt={'For smoothing pdgm into null continuum:'};
         def={['[' int2str(daniel2) ']']};
         tit1='Input the Daniel Filter Lengths';
         lineNo=1;
         answer=inputdlg(prompt,tit1,lineNo,def);
         daniel2= str2num(answer{1});
      end
      wnull = danwgtn(daniel2);
      mm=length(wnull); % span of daniel filter
      igo = n+1+(mm-1)/2;
      istop=2*n+(mm-1)/2;
      ynull=conv(yext,wnull); % smoothed periodogram == null line
      hp7b=plot(ff,y,'.g',ff,ysm,'-b',ff,ynull(igo:istop),'-m');
      set(hp7b(1),'Color',cord(3,:));
      set(hp7b(2),'Color',cord(1,:));
      set(hp7b(3),'Color',cord(2,:));
      
       % Plot bandwidth bar.  Change next 2 lines as needed for formatting
      xlim7=get(gca,'XLim'); 
      ylim7=get(gca,'YLim');
      delx = diff(xlim7);
      dely = diff(ylim7);
      ybwgo = ylim7(1)+0.83*dely;  % bar 8/10 of way to top of figure
      yinc = dely/150; % vertical width of tick at ends ob bw bar
      xbwgo=xlim7(1)+delx/2 - bw/2;
      xbwsp = xbwgo+bw;
      hlineb1=line([xbwgo xbwsp],[ybwgo ybwgo]); % the bw bar itself (horiz line)
      hlineb2=line([xbwgo xbwgo],[ybwgo-yinc/2 ybwgo+yinc/2]); % left side vert line
      hlineb3=line([xbwsp xbwsp],[ybwgo-yinc/2 ybwgo+yinc/2]); % right side vert line
      set(hlineb1,'Color',1-get(gca,'Color'));
      set(hlineb2,'Color',1-get(gca,'Color'));
      set(hlineb3,'Color',1-get(gca,'Color'));
      
      text(xbwgo+(xbwsp-xbwgo)/2,ybwgo+yinc/2,'BW',...
         'VerticalAlignment','bottom','HorizontalAlignment','center');
            
      xlabel('Frequency (yr^{-1})');
      ylabel('Relative Variance');
      grid;
      set(gca,'Position',[0.1300    0.1100    0.7750    0.78]);
      title({['Smoothed Pdgm of ' name ';  ' stryr],...
            ['Spans: [' int2str(daniel1) '], [' int2str(daniel2) ']']});
      set(hp7b(2),'LineWidth',1.5);
      set(hp7b(3),'LineWidth',1.5);
   elseif kmen1==2;
      kwh1=0;
      ynull = ynull(igo:istop);
   end
end



%********* CONDIDENCE INTERVAL AROUND SPECTRAL ESTIMATES  ******
%
% A chi-sq test is used.  Degrees of freedom depends on (1) original 
% series length,  (2) padded length, (3) proportion of data tapered,
% and (4) filter weights.

% Compute the factor gsq, that depends on the above four items.  First
% handle tapering.  See Bloomfield, p. 194.
ufact=0.5*(128-93*mtap1)/((8-5*mtap1)*(8-5*mtap1));

% Variance depends on filter weights.
% Compute sum of squares of weights.
wfact=w'*w;

% Bring adjustment for padding (Bloomfield, p. 192) into computation of
% complete variance adjustment for spectral estimates (Bloomfield, p.
% 195, eq. 12).
gsq = ufact*wfact*padlen/nyra;

% Compute degrees of freedom for spectral estimates
vdf = 2/gsq;   



%********************* FINISHED SPECTRAL PLOT

% find the row number of freq in ff closest to key freq
if strcmp(khypo,'Yes');
   d = ff - keyfreq;
   [ymin,imin]=min(abs(d));
else
   ymin=[]; imin=[];
end



kwhl2=1;
while kwhl2==1;
   kmen2=menu('Choose 1',...
      'Make or revise final spectral plot',...
      'Satisfied -- quit program');
   
   if kmen2==1; % make or revise
      kmen3=menu('Choose one',...
         'Spectrum with 95% CI',...
         'Spectrum with 99% CI',...
         'Spectrum without CI');
      if kmen3==1;
         chi1 = chi2inv(0.975,round(vdf)); 
         chi2 = chi2inv(0.025,round(vdf));
         strci='95% CI';
         % Construct confidence intervals around the spectral estimates  using
         % the chi-sq distribution (Bloomfield p. 196)
                  
      elseif kmen3==2;
         chi1 = chi2inv(0.995,round(vdf));
         chi2 = chi2inv(0.005,round(vdf));
         strci='99% CI';
      elseif kmen3==3;
         chi1=[];
         chi2=[];
         strci=[];
      end
      if kmen3==1 | kmen3==2;
         figure(8);
         ysmlo = vdf/chi1 * ysm;
         ysmhi = vdf/chi2 * ysm;
         hp8=plot(ff,ysm,'-b',ff,ynull,'-m',ff,ysmlo,'-m',ff,ysmhi,'-m');
         set(hp8(1),'Color',cord(1,:));
         set(hp8(2),'Color',cord(2,:));
         set(hp8(3),'Color',cord(2,:));
         set(hp8(4),'Color',cord(2,:));
                     
          % Plot bandwidth bar.  Change next 2 lines as needed for formatting
          xlim8=get(gca,'XLim'); 
          ylim8=get(gca,'YLim');
          delx = diff(xlim8);
          dely = diff(ylim8);
          ybwgo = ylim8(1)+0.83*dely;  % bar 8/10 of way to top of figure
          yinc = dely/150; % vertical width of tick at ends ob bw bar
          xbwgo=xlim8(1)+delx/2 - bw/2;
          xbwsp = xbwgo+bw;
          hb1=line([xbwgo xbwsp],[ybwgo ybwgo]); % the bw bar itself (horiz line)
          hb2=line([xbwgo xbwgo],[ybwgo-yinc/2 ybwgo+yinc/2]); % left side vert line
          hb3=line([xbwsp xbwsp],[ybwgo-yinc/2 ybwgo+yinc/2]); % right side vert line
          set(hb1,'Color',1-get(gca,'Color'));
          set(hb2,'Color',1-get(gca,'Color'));
          set(hb3,'Color',1-get(gca,'Color'));
          
          text(xbwgo+(xbwsp-xbwgo)/2,ybwgo+yinc/2,'BW',...
             'VerticalAlignment','bottom','HorizontalAlignment','center');
          
          % Complete figure
          xlabel('Frequency (yr^{-1})');
          ylabel('Relative Variance');
          set(gca,'Position',[0.1300    0.1100    0.7750    0.78]);
          title({['Smoothed Pdgm of ' name ';  ' stryr],...
                ['Spans: [' int2str(daniel1) '], [' int2str(daniel2) ']']});
          set(hp8(1),'LineWidth',1.5);
          set(hp8(2),'LineWidth',1.5);
          legend('Spectrum','Null Continuum',strci);
          
          % Add vertical line at test frequency
          if strcmp(khypo,'Yes');
             hlinetest=line([keyfreq keyfreq],ylim8);
             set(hlinetest,'Color',cord(3,:));
          end
          

      else
         figure(8);
         hp8=plot(ff,ysm,'-b', ff,ynull,'-m');
         set(hp8(1),'Color',cord(1,:));
         set(hp8(2),'Color',cord(2,:));
         % Plot bandwidth bar.  Change next 2 lines as needed for formatting
          xlim8=get(gca,'XLim'); 
          ylim8=get(gca,'YLim');
          delx = diff(xlim8);
          dely = diff(ylim8);
          ybwgo = ylim8(1)+0.83*dely;  % bar 8/10 of way to top of figure
          yinc = dely/150; % vertical width of tick at ends ob bw bar
          xbwgo=xlim8(1)+delx/2 - bw/2;
          xbwsp = xbwgo+bw;
          hb1=line([xbwgo xbwsp],[ybwgo ybwgo]); % the bw bar itself (horiz line)
          hb2=line([xbwgo xbwgo],[ybwgo-yinc/2 ybwgo+yinc/2]); % left side vert line
          hb3=line([xbwsp xbwsp],[ybwgo-yinc/2 ybwgo+yinc/2]); % right side vert line
          set(hb1,'Color',get(gca,'Color'));
          set(hb2,'Color',get(gca,'Color'));
          set(hb3,'Color',get(gca,'Color'));
          
          text(xbwgo+(xbwsp-xbwgo)/2,ybwgo+yinc/2,'BW',...
             'VerticalAlignment','bottom','HorizontalAlignment','center');
          
          % Complete figure
          xlabel('Frequency (yr^{-1})');
          ylabel('Relative Variance');
          set(gca,'Position',[0.1300    0.1100    0.7750    0.78]);
          title({['Smoothed Pdgm of ' name ';  ' stryr],...
                ['Spans: [' int2str(daniel1) '], [' int2str(daniel2) ']']});
          set(hp8(1),'LineWidth',1.5);
          set(hp8(2),'LineWidth',1.5);
          
        
      end
      
   elseif kmen2==2;
      kwhl2=0;
   end
end

%************** ascii output

kascii = questdlg('Ascii summary file?')
if strcmp(kascii,'No') | strcmp(kascii,'Cancell');
   % no action needed
else
   [file1,path1]=uiputfile('info*.dat','Ascii info to go here');
   pf1=[path1 file1];
   
   head1=' N  Frequency   Period (yr)  Spectrum';
   
   fid1=fopen(pf1,'w');
   fprintf(fid1,'%s\n',name);
   fprintf(fid1,'%s\n',blanks(5));
   fprintf(fid1,'Analysis Period: %4.0f-%4.0f\n',yrsa);
   fprintf(fid1,'Time Series Length = %4.0f yr\n',nyra);
   fprintf(fid1,'Padded Length = %4.0f yr\n',padlen);
   fprintf(fid1,'Frequency spacing = %8.6f \n',1/padlen);
   fprintf(fid1,'%s\n',blanks(5));
   
   fprintf(fid1,'Daniell spans for null line: %s\n',int2str(daniel2));
   fprintf(fid1,'Daniell spans for spectrum: %s\n',int2str(daniel1));
   fprintf(fid1,'Bandwidth for Spectral Estimates = %8.6f \n',bw);

fprintf(fid1,'%s\n',blanks(5));

if strcmp(khypo,'Yes');
   fprintf(fid1,'Frequency of hypothesized periodicity = %8.5f\n',keyfreq);
   fprintf(fid1,'   Frequency for test: %8.5f\n',ff(imin));
   fprintf(fid1,'Period of hypothesized periodicity = %8.5f yr\n',keypd);
   fprintf(fid1,'   Period for test: %8.5f\n',1/ff(imin));
   fprintf(fid1,'%s\n',blanks(5));

else
end


  
   fprintf(fid1,'%s\n',strv1);
   fprintf(fid1,'%s\n',strv2);
   fprintf(fid1,'%s\n',strvpct);
   
   fprintf(fid1,'%s\n',blanks(5));

   if ~isempty(strci);
      fmt1 = '%3.0f  %6.5f  %8.2f       %10g\n';
      fprintf(fid1,'Peaks significant at %s\n',strci);
      Lsig = ynull<ysmlo;
      sum1 = sum(Lsig);
      if sum1>0;
         ff = ff(Lsig);
         ysm=ysm(Lsig);
         ysmlo=ysmlo(Lsig);
         per = 1 ./ ff;
         
         fprintf(fid1,'%s\n',head1);
         fprintf(fid1,'%s\n',blanks(5));
         
         for nn = 1:sum1;
            fprintf(fid1,fmt1,nn,ff(nn),per(nn),ysm(nn));
         end
      end
   end
end


fclose all
         
         
   
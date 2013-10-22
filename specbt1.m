function G=specbt1(z,M,f,prevwind,sername,kopt)
% specbt1:  spectrum of a single time series  by Blackman Tukey method
% G=specbt1(z,M,f,prevwind,sername,kopt);
% Last revised 9-16-99
%
% Compute and optionally plot the sample spectrum with 95% confidence bands. 
% You can interactively control the spectral smoothing through M, the length 
% of the smoothing window.  Figure includes white noise spectral line and 
% bandwidth bar.
%
%*** IN 
%
% z (mz x 1)r  the time series 
% M (1 x 1)i  length of lag window used in calculations.  
% f (mf x 1)r  frequencies (cycles/yr, in range 0-0.5) for spectral estimates
% prevwind(1 x 1)i  previous figure window % sername(1 x ?)s   name of series:  for annotating title
% kopt (1 x 1)i options
%    kopt(1): whether or not to plot within function 
%      ==1 interactive plotting (linear y) to select final value of M
%      ==2 no plotting; just compute spectrum given initial M
%      ==3 interactive plotting using semilogy plot
%
%*** OUT 
%
% G (mG x 5)r   spectrum and related quantities. By column:
%   1- frequency
%   2- period (yr)
%   3- spectral estimate
%   4- lower 95% confi limit
%   5- upper 95% conf limit
%
% Graphics output:  
% Figure 1.  Spectrum, with 95% confidence interval.  Horizontal line marks
% spectrum of white noise with same variance as z. Horizontal bar marks
% bandwidth of Tukey window used to smooth autocovariance function
%
%*** REFERENCES 
% 
% Chatfield, C. 1975.  Time series analysis, theory and practice.  
% Chapman and Hall, London, 263 p.  
%
%*** UW FUNCTIONS CALLED
%
% specbt1/calcspec -- subfunction that does most of the meaty calculations
%
%*** TOOLBOXES NEEDED
%
% statistics
% system identification
%
%*** NOTES
%
% More info on the Input:
%
%   M:  Suggest some value on the order of 1/3 to 1/20 of the time series length; 
%   increased M leads to more detail in spectrum, but greater variance of estimates
%
%   f:  Example: [0:.01:.5] gives estimates at freq 0, 0.01, 0.02, ...
%
%   prevwind:  included to allow you to avoid overwriting windows you might have 
%   in the calling program.  If you want spectrum in window 1, just set 
%   prevwind to 0.
%
% Spectrum computed by MATLAB function spa.m, which implements Blackman-Tukey 
% method for spectral estimate.
%
% 95% Confidence interval computed by equation in Chatfield (1975, p 150).  
% Degrees of freedom, following Chatfield, computed as 2.67N/M, where N
% is the sample size and M is the length of the lag window.  Critical points
% of Chi-squared distribution (.025 and .975) computed by call to MATLAB
% function chi2inv.m.
%
% Bandwidth defined in Chatfield (p. 154).  Bandwidth for Tukey window is
% (8 * pi)/(3 * M) in radians.  Using conversion omega=2*pi*f, where omega
% is frequency in radians and f in in cycles per time unit, bandwidth in 
% cycle/unit time is  bw=4/(3M).
%
% Note that spa is in the 'system identification toolbox'. The function desription
% says a Hamming window is used to smooth the covariance function.  The Hamming
% window as defined in Ljung (1987, p. 155) is the same as what Chatfield(1975,
% p. 140) calls the Tukey window.


nyr = length(z);
defM = min([30 round(nyr)/10]); % default length of lag window

% Subtract mean
z = z - mean(z);

% Initial computation and plotting of spectrum
G=calcspec(z,M,f,prevwind,sername,kopt);

if kopt(1)==2;
   return
end


kinit=1;  % initial setting for default bandwidth

k2=1;
k1 = 1;
while k1~=0;
   k2=menu('Choose 1',...
      'Change Length of Lag Window',...
      'Change y-axis scale',...
      'Accept Spectrum and Quit');
   
   
   switch k2
   case 1
      prompt={'Enter the value for M:'};
      def={num2str(M)};
      titl='Change Lag Window';
      lineNo=1;
      M=inputdlg(prompt,titl,lineNo,def);
      M=str2num(M{1});
      
      G=calcspec(z,M,f,prevwind,sername,kopt);
   case 2; % change scale
      if kopt(1)==1;
         kopt(1)=3;
      elseif kopt(1)==3;
         kopt(1)=1;
      end;
      G=calcspec(z,M,f,prevwind,sername,kopt);
     
            
   case 3
      k2=0;
      k1=0;
   otherwise
      
   end
   
end

function G=calcspec(z,M,f,prevwind,sername,kopt)
% Convert frequency
N=length(z);
w = f * 2 * pi;  
% Compute spectrum
g = spa(z,M,w); % freq in col 1, spectrum in col 2 , st dev in col 3
g(1,:)=[];  % remove first row, which holds variable identifiers

% Frequency from radian units to cycles per year
f = g(:,1)/(2*pi);

% Compute degrees of freedom for Tukey Window
df = round(2.67*N/M);  % see Chatfield (1975, p. 150)

% Find critical points of Chi Squared for 95% confidence interval
c = chi2inv([.025 .975],df);

% Compute confidence interval
clow = df * g(:,2)/c(1);
chi = df * g(:,2)/c(2);

% Store spectal info
nf = length(f);
f1 = f(2:nf);
p1 = 1 ./ f1;
p2 = [inf; p1];
G=[f p2 g(:,2) clow chi];


if kopt(1)==2;
   return
end


% Plot spectrum
figure(prevwind+1);
if kopt(1)==1;
   hspec=plot(f,g(:,2),...
      f,clow,'-r',...
      f,chi,'-r');
elseif kopt(1)==3;
   hspec=semilogy(f,g(:,2),...
      f,clow,'-r',...
      f,chi,'-r');
end;
set(hspec(2),'Color',[1 .5 1]);
set(hspec(3),'Color',[1 .5 1]);
xlabel('Frequency (yr^{-1})');
ylabel('Relative Variance');
legend('Spectrum','95% Conf. Interval');

% White noise line
hline=line([min(f) max(f)],[var(z) var(z)]);
set(hline,'linestyle',':');

 %Compute and plot bandwidth
bw = 4/(3*M);  % see Ljung (1987, p. 155), Chatfield (1975, p. 154)
% Note that the "Hamming" window as defined in system ident toolbox and Ljung
% is identical to Chatfield's Tukey Window
% Put bandwidth bar 0.3 from top of figure, centered
ylims = get(gca,'Ylim');
yrng = abs(ylims(2)-ylims(1));
if kopt(1)==1;
   ypoint = ylims(2)-0.3*yrng;
   line([0.25-bw/2 0.25+bw/2], [ypoint ypoint]);
   line([0.25-bw/2  0.25-bw/2],[ypoint+yrng/100 ypoint-yrng/100]);
   line([0.25+bw/2  0.25+bw/2],[ypoint+yrng/100 ypoint-yrng/100]);
   htt = text(0.25,ypoint+yrng/100,'Bandwidth');
elseif kopt(1)==3;
   dlog=log10(ylims(2))-log10(ylims(1));
   dlog5=dlog/5;
   dfact = 10^(-dlog5);
   ypoint = dfact*ylims(2);
   line([0.25-bw/2 0.25+bw/2], [ypoint ypoint]);
   htt = text(0.25,ypoint,'Bandwidth');


end;


set(htt,'HorizontalAlignment','Center','VerticalAlignment','bottom');

% Build title
str1 = sprintf('%4.0f',N); % Length of time series
str2 = sprintf('%4.0f',M); % Width of lag window
str3 = sername; 
title(['Spectrum of ' sername ';  Lag Window \itM/N=\rm 'int2str(M) '/' int2str(N)]);


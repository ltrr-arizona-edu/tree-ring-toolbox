function [tcy,f,S,M,m1,m2]=evolspec(x,yrgo,tshift,tslab,padlen,bw,figwin)
%
% Evolutionary spectral analysis 
%
% D. Meko 5-30-95
%
%********* IN ARGS ***************
%
% x (mx x 1)r	time series
% yrgo (1 x 1)i	start year of x
% tshift( 1 x 1)i  number of years  to offset slabs
% tslab(1 x 1)i width (number of years) of time slabs
% padlen (1 x 1)i padded length of time slabs -- must be a power of 2
%    and must be at least as large as tslab
% bw 1 x 1)r   bandwidth of the spectral window, must be 
%	between 0 and 0.5
% figwin (1 x 1)i  first figure window to use (will use this and next higher)
%
%****************  OUT ARG ***********
%
% tcy (col vector)i  central years of the time slabs
% f(padlen/2 x 1)r   frequency points for plotting spectra
% S (mS x nS)r	estimated spectra;  each col is a spectra for a time
%	slab;  each row is a central frequency of spectral estimate
% M (1 x 1)i the "M" parameter value used by etfe.m
% m1(1 x 1)i  length of the Hamming filter used to smooth pdgm
% m2(1 x 1)r  length of the rectangular filter that would give
%	the desired bandwidth bw
%
%***************** USER-WRITTEN FUNCTIONS CALLED ** none
%
 
%*********** METHOD ****************************************
%
% Spectra are estimated for offset segments (slabs) of the time series.
% Length and offset of slabs are controlled with tslab and tshift.
% Resolution is controlled with the bandwidth (bw) and the padding
%	length (padlen)
% The matrix (S) of the spectra is plotted with a color map and a
%	contour plot.  The abscissa is the central year of the slabs.
%	the ordinate is frequency (cycles/year), going from 0 to 0.5.
%	A bandwidth "bar" is placed with the mouse wherever you want
%	on the plot once the plot is finished on the screen.
%
% Details of the estimation of spectra for individual slabs:
%   
%	The mean is first subtracted.
%	The series is then padded to specified length padlen
%	The FFT and periodogram is computed
%	The periodogram is smoothed with a Hamming frequency window
%		with the desired bandwidth (see below)
%
% Bandwidth of frequency window;  selection of Hamming frequency window
%
%   Following Chatfield (p. 154), the bandwidth is defined as the width
%   of the "ideal" rectangular spectral window which would give an 
%   estimator with the same variance as the selected Hamming window.
%   The desired bandwidth (bw) is specified as an input argument.  
%   The equivalent width, m2, of the rectangular filter is bw*padlen, where
%   padlen is the padded length of the time series. 
%
%   Hamming filters with increasing width, beginning with m1=m2, are
%   computed using the built-in function hamming.m.  The filter weights
%   are scaled so that they sum to 1.0.  The sum-of-squares for each
%   computed Hamming filter are computed, and the filter whose SOS
%   is closest to that of the rectangular filter with width m2 is
%   identified.  This is initial selected width of the Hamming filter 
%   to be used.  
%
%   This initial m1 is adjusted slightly, if necessary, to be
%   compatible with the "M" argument actually used by matlab function
%   etfe.m to compute the spectrum.  As noted in the SI toolbox, etfe
%   smooths the periodogram with a Hamming filter of length 
%   m1 = fix(padlen/M) + 1
%   Thus the Hamming-filter length m1 is not directly passed as an 
%   input to etfe.m, but is computed within etfe.m from the 
%   parameter "M".  Evolspec.m handles the problem by
%   using the starting value of m1 described in the pgh above, and then
%   computing M from 
%      M = fix(padlen/(m1-1))
%   The final value of m1 is then computed as 
%   m1 = fix(padlen/M) + 1
%
% References
%
% Ljung, Lennart, 1987.  System Identification: 
% Theory for the User.  Prentic-Hall, Inc., 519 pp.  
%  Section 6.4 "Spectral Analysis" (p. 151) describes smoothing the
%  empirical transfer-function estimate (etfe) as a method of
%  getting the spectrum.  The Hamming frequency window is discussed
%  on p. 155.  
%
% Ljung, Lennart, 1991. System Identification Toolbox, User's Guide.
% The Mathworks, Inc.  
%  Evolspec.m uses the Matlab function etfe, which is described on 
%  p. 2.33.  The statement that the length of the Hamming window is equal
%  to "the number of data points in z divided by M, plus one" is
%  incorrect, except for the case when the number of data points in 
%  z happens to equal the specified length of the padded series.
%
% Chatfield, C., 1975.  The Analysis of Time Series: Theory and 
% Practice.  John Wiley & Son, 263 pp.  
%  Concept of bandwidth is discussed on p. 154.  Tukey window, 
%  defined by eqn on p. 140, is identical to Hamming window, as
%  defined by Ljung (1987, p. 155).
%
% Bloomfield, Peter, 1976.  Fourier Analysis of Time Series: an 
% Introduction.  John Wiley & Sons, 258 pp.  
%  The raw periodogram as a method for estimating
%  the spectrum.  Means and variances of smoothed spectra, and 
%  necessary adjustments for padding (p. 224).  Proportionality of
%  variance of the spectral estimates on the sum-of-square of
%  filter coefficients (p. 224).  This proportionality is used in 
%  evolspec.m to find the length of the Hamming filter for 
%  smoothing the periodogram such that the desired bandwidth of
%  spectral estimates is reached.
%  


%***************** ADJUSTING COLOR AND CONTOUR PLOTS
%
% Grid lines might get in the way in the colormap. Use
% "shading" command to vary how lines show
%   shading faceted, interp, flat
%       faceted is default, and gives black grid lines
%	flat gets rid of lines, but image can be jagged
%	interp gives a nice smooth image with no grid lines,
%		but takes longer
%
% See commented out sections at end of this function for 
%   display variations that can be turned back on.  For example,
%   contour lines can be labeled on the contour plot. A neat
%   super-looking color plot results from overlaying contours
%   after getting rid of grid on the color map -- again, see
%   commented out section.
%
% Zooming can be done with the Matlab zoom command in the 
%   figure window.  This is a way to get a plot of a restricted
%   frequency or time range of the full evolutionary spectral plot
%
%
%****************  PLOTTING RESULTS WITH OUTSIDE SOFTWARE  ********
%
% This is no problem, since the output arguments (see above) include
% everything needed to produce the colormap and contour plot.  
%
%
%***************** TIPS ON RUNNING  ***********************************
%
% User will need to point with mouse to point on plots at which to
% position top of "bandwidth" bar. No screen prompt appears -- you
% just need to know this.
%
% After bar is plotted on figure, press Return.  When both figs 1 and 2
% are complete, you can toggle back through the windows to look em over.
%
% Shortcut to repeatedly typing in function call is to store a  command
% string as one of the variables with the data in the in the .mat file.
% For example, make a string variable "doit" as follows:
%  doit= '[tcy,f,S,M,m1,m2]=evolspec(x,yrgo,tshift,tslab,padlen,bw);'
%
% Save doit along with other variables in the .mat file.  Then, after
% loading the .mat file, simply key in:
%  eval(doit)


% Size and check arguments
[mx,nx]=size(x);
if ~ (mx>=100 & nx==1),
   error('x must be col vector of min length 100')
end
if tslab >= mx,
   error('Time slab must be shorter than length of x')
end
if rem(tslab,2)==0
   error('Time slab must have odd number of years')
end
if tshift <=0
	error('Time shift must be at least 1 year')
end
if bw>=0.5,
  error('Required: bw<0.50')
end


% Check that padlen is a power of 2 and is at least as large as tslab
if padlen<tslab, error('Required: padlen>= tslab'), end;
np2 = nextpow2(padlen);
if 2^np2 ~= padlen, 
   error('padlen must be power of 2')
end 


% Compute some miscellaneous quantities
yrsp = yrgo + mx -1;  % ending year of time series x
tey = tslab+yrgo-1:tshift:yrsp;   % ending years of slabs
nslabs = length(tey); % number of slabs
tcy = (ceil(tslab/2)+yrgo-1):tshift:yrsp;  % central years of slabs
tcy = tcy(1:nslabs); % truncate so that all tcy are for valid
	% segments (i.e., do not over-reach ends of time series)


% Compute row-index matrix into x
t1 = 1:tshift:mx;  % row index to beginning years of slabs
t1=t1(1:nslabs);
t2 = t1(ones(tslab,1),:); % expand t1 to matrix
j1 = (0:tslab-1)';  % col vector
j2 = j1(:,ones(nslabs,1)); % expand j1 to matrix
T = t2+j2;  % Index to rows of x for pulling out slabs


nfft=padlen;  % will be computing the padlen-length fft
N = padlen/2;  % number of frequency points for plotted spectra

% Compute length of rect filter with desired bandwidth
m2 = padlen* bw;  % note that the sum-of-squares of this rect filter is 1/m2, and
   % that m2 is 1 divided by the sum of squares of the weights
m2 = fix(m2);  % nearest integer (rounded towards zero)



%***************** Find Hamming filter that will give estimate with same variance
% as that rectangular window
jvect = m2:1:min([N 2*m2]);  % will try windows of these lengths
njvect = length(jvect); % will try this many window lengths
M2 = repmat(NaN,njvect,1);
J = repmat(NaN,njvect,1);
for j = 1:njvect;
   jj=jvect(j); % length of the window
   ha = hamming(jj); % hamming filter
   ha = ha/sum(ha);  % hamming filter with weights summing to 1.0
   M2(j) = 1/ sum(ha .* ha);  % 1 divided by the sum of squares of the hamming weights
   J(j)=jj; % length of the Hamming filter
end
% The desired Hamming filter will be the one with the sos of filter weights closest to
% 1/m2, or equivalently, with 1/sos(filter weights) closest to m2
d1 = abs(M2-m2);
[d2,I]=min(d1);
m1 = J(I);
if(m1<3),
  error('Bandwidth too narrow for this padded length')
end


% Compute the "lag-window" length M that corresponds to the 
% frequency window length m1.  See etfe.m description in 
% syst ident toolbox
M = fix(padlen/(m1-1));
m1 = fix(padlen/M)+1;  % computed length of Hamming freq window
	% this is the length of window etfe.m will actually smooth
	% periodogram with

% Get Hamming filter for frequency window
ha = hamming (m1);

% Run a dummy spectral analysis to get the plotting-point frequencies
% that will apply to all the spectra
gdum = etfe(randn(tslab,1),M,N);
gdum(1,:)=[];
f = gdum(:,1)/(2*pi); % abscissa will run from 0 to 0.5 in units 1/year

c = NaN;
S = c(ones(N,1),ones(nslabs,1));


% Compute spectra
for n = 1:nslabs
   t = T(:,n);  % row index in x for this slab
   z = x(t); 
   z = z - mean(z);  % convert slab to zero-mean
   G = etfe(z,M,N); % estimate spectra
   G(1,:)=[];
   g = G(:,2);
   S(:,n)=g;
end


% Figure figwin will have the color map
figure(figwin);
h1=pcolor(tcy,f',S);
[xx1 yy1]=ginput(1);
yyy1=yy1-bw;
hl1 = line([xx1 xx1],[yy1 yyy1]);
htxt1=text(xx1,yy1+0.03,'BW','HorizontalAlignment','center','color',[0 0 0],...
   'FontSize',10);
set(hl1,'Linewidth',3,'color',[0 0 0])
title(['EVOLUTIONARY SPECTRUM, WINDOW= ',int2str(tslab),...
	' YR,  Shift = ',int2str(tshift),' YR'])
xlabel('Year')
ylabel('Frequency (yr^{-1})');
hylab = get(gca,'ylabel');
set(hylab,'erasemode','none');
shading interp
%hold on
%contour(tcy,f',S,10,'k')
%hold off
pause


figure(figwin+1);
h2=contour(tcy,f',S);
[xx1 yy1]=ginput(1);
yyy1=yy1-bw;
hl1 = line([xx1 xx1],[yy1 yyy1]);
text(xx1,yy1+0.03,'BW','HorizontalAlignment','center','Color',[0 0 0]);
set(hl1,'Linewidth',3,'color',[0 0 0])
title(['EVOLUTIONARY SPECTRUM, WINDOW= ',int2str(tslab),...
	' YR,  Shift = ',int2str(tshift),' YR'])
xlabel('Year')
ylabel('Frequency (yr^{-1})');

figure(figwin);




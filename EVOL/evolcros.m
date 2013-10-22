function [tcy,f,C,P,bwout,p,spans,w,df]=evolcros(x,y,yrgo,tdat,bw,ns,taper)
% usage: [tcy,f,C,P,bwout,p,w]=evolcros(x,y,yrgo,tdat,bw,ns,taper)
%
% Evolutionary cross-spectral analysis 
%
% D. Meko 7-27-95
%
%********* IN ARGS ***************
%
% x (mx x 1)r	time series # 1
% y (my x 1)r	time series # 2
% yrgo (1 x 1)i starting year of time series x and y
% tdat (1 x 3) time-slice information:
% 	tshift --   number of years  to offset slabs
% 	tslab --  width (number of years) of time slabs
% 	padlen -- padded length of time slabs -- must be a power of 2
%    and must be at least as large as tslab
% bw (1 x 1)r   desired bandwidth of the spectral window, must be 
%	between 0 and 0.5
% ns (1 x 1)i  number of passes of Daniell filter desired in smoothing
% taper (1 x 1)r total proportion of segments of series to be tapered
%
%****************  OUT ARG ***********
%
% tcy (col vector)i  central years of the time slabs
% f(padlen/2 x 1)r   frequency points for plotting spectra
% C (mC x nC)r	estimated squared coherency  each col is a spectra for a time
%	slab;  each row is a central frequency of spectral estimate
% P (mP x nP)r estimated phase -- organized like C
% bwout (1 x 1)r bandwidth of final convoluted Daniell filter
% p (1 x 3)r  90%, 95% and 99% confidence levels for C
% spans (1 x ns) % spans of individual Daniel filters
% w (mw x 1)r  weights of final Daniell filter
% df (1 x 1)r degrees of freedom for spectral estimates
%

%*********** METHOD ****************************************
%
% Computation of spectra and cross-spectra follow the book 
% "Fourier Analysis of Time Series:  An Introduction", 
% 1976, Peter Bloomfield, 258 pp, John Wiley & Sons
%
% Bloomfield includes listings of Fortran programs for computing
% spectra and cross-spectra. Equations, but not Fortran code are
% included for confidence bands on squared coherency.
%
% Extensions beyond Bloomfield include the following:
%	1.  Automatic selection of spans of Daniell filter to give
%		spectral estimates with desired bandwidth. The user just
%		needs to specify the bandwidth and the number of passes of
%		the filter.  The final filter is computed automatically.
%  2.	More flexible padding with zeros option.  
%		Series can be padded with zeros to a length greater than
%		the first power of two above the series length.  This extension
%		might be useful in getting spectral estimates centered on
%		frequencies of a-priori interest. 
%  3.	The "evolutionary" part of the problem.  This amounts to
%		repeatedly calling for cross-spectral computations on different
%		segment of the time series.
%
%
%******************** SUMMARY OF STEPS ******************************
%
% Size and check the input arguments
%
% Compute the Daniel filter, its sum-of-squares, and other properties
%
% From the filter weights, compute the theoretical 95% and 99% 
% confidence bands that will apply to the squared-coherency estimates
%
% Set up the index matrices to get segments from the time series
%
% Loop over segments
%
%		Get segment of data for each series x and y
%		Subtract segment mean; taper and pad, depending on options
%		Compute Fourier transforms and raw periodograms
%		Smooth periodograms to get estimated spectra
%		Compute cross-periodgram from Fourier transforms
%		Smooth Cross periodogram and use it and the smoothed
%			spectra to get estimate of cross-spectrum
%		Compute squared coherency from amplitude of cross-spectrum
%		Store squared coherency and phase for segment
%
% Make color map and contour map of squared coherency and phase
%
%******************* USER - WRITTEN FUNCTIONS NEEDED ****************
%
% danbw -- get daniell filter from number of spans and bandwidth
% danwgtn -- get daniell filter convoluted n times
% cohprob1 -- .90,.95, .99 prob points for dist of squared coherency
% crspec1 -- compute estimated cross-spectrum
%
%****************  NOTES ********************************************
%
%
% Squared coherency and phase have been checked against values 
% output by SPLUS function spec.pgram on same data sets
% 
% Spectra are estimated for offset segments (slabs) of the time series.
% Length and offset of slabs are controlled with tslab and tshift.
% Resolution is controlled with the bandwidth (bw) 
%
% The matrices C and P hold squared coherency and phase for the segments
%	(slabs) of the series.  These matrices are displayed with a color map
%	and a contour plot.  The abscissa is the central year of the slabs.
%	The ordinate is frequency (cycles/year), from 0 (infinite period) to
%	to 0.5 (Nyquist period of 2 years).	A bandwidth "bar" is placed with 
%	the mouse wherever you want on the plot once the plot is finished on the 
%	screen.
%
% Bandwidth of frequency window;  selection of the appropriate Daniell 
% filter to smooth the periodograms and cross-periodogram
%
%   Following Chatfield (p. 154), the bandwidth is defined as the width
%   of the "ideal" rectangular spectral window which would give an 
%   estimator with the same variance as the selected spectral window.
%   The desired bandwidth (bw) is specified as an input argument. Also
%	 specified is the number of passes, ns, for the Daniell filter. The 
%	 final filter is computed by convoluting a basic Daniell filter of
%	 length m with itself, then convoluting the result with the basic
%	 filter, etc., each time producing a filter with a more bell-shaped
%	 appearance.  If m is the length of the elemental Daniell filter,
%	 one convolution gives a filter with length m+(m-1).  Two convolutions
%	 give a filter with length m + (m-1)+(m-1).  N covolutions yield a
%	 filter length m+N(m-1).  N convolutions is equivalent to a filter
%	 with ns = N+1 spans.  Thus, the final filter length is
%	 L = m+ns(m-1).
%
%	 Given the bandwidth bw and the padded length of the series, one
%	 can compute the length of the rectangular filter needed.  The Daniell
%	 filter length L must be greater than this. The program begins by
%	 setting L equal to the length of the ideal rectangular filter, 
%	 computing an initial Daniell length m from
%		L = m+ns(m-1),
%	 computing the convoluted Daniell filter weights for this m and ns,
%	 checking the sum of squares of these weights against the sum of 
%	 squares of the ideal rectangular filter, and then iteratively 
%	 increasing m by 2 until the sum-of-squares is closest to the that
%	 of the ideal rectangular filter.
%
% Additional Reference
%
%
% Chatfield, C., 1975.  The Analysis of Time Series: Theory and 
% Practice.  John Wiley & Son, 263 pp.  
%  Concept of bandwidth is discussed on p. 154.  Tukey window, 
%  defined by eqn on p. 140, is identical to Hamming window, as
%  defined by Ljung (1987, p. 155).
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



% Save doit along with other variables in the .mat file.  Then, after
% loading the .mat file, simply key in:
%  eval(doit)


t0=clock; % to check elapsed time until begin plotting
% Can put dbstop at line 341 to check time1


close all


% Size and check time series
[mx,nx]=size(x);
[my,ny]=size(y);
if (~ (mx>=100 & nx==1)) & (~(my>=100 & ny==1)),
   error('x and y must be col vectors of min length 100')
end

% Rename variables for convenience
tshift=tdat(1);
tslab=tdat(2);
padlen=tdat(3);

% Check for valid settings on time slabs and shift
if tslab >= mx,
   error('Time slab must be shorter than length of x')
end
if rem(tslab,2)==0
   error('Time slab must have odd number of years')
end
if tshift <=0
	error('Time shift must be at least 1 year')
end
if bw<=0 | bw>=0.5,
  error('Required: bandwidth must be GT 0 and LT 0.5')
end

% Check that taper proportion is reasonable. Set maximum for now as
% 0.2 (20 percent) of the total slab length.  This means no more
% than 10 percent on each end.
if ~(taper >=0 & taper <=0.2),
	error('taper should be in range 0 to 0.2')
end 


% Check that padlen is valid.  Padded length must be a power of two,
% must be at least as long as the time slabs, and must
% be no more than 32 times the time slabs.  The "32" has no special
% meaning. Who in their right mind would want to pad a 1000 year series
% to length greater than 32*1024 anyway?
if padlen<tslab, error('padlen shorter than time slabs'), end
ncheck = nextpow2(tslab);
next2 = 2 ^ncheck;
if padlen>32*next2, 
	error('padlen more than 32X next power of two greater than tslab')
end
n2 = next2*[1 2 4 8 16 32];
if ~any(n2==padlen)
	error('padlen not a power of 2')
end


% Check that bandwidth is reasonable
if bw <=0 | bw>=0.5,
	error('bandwidth (bw) must be greater than 0 and less than 0.5')
end


% Check that number of passes of the Daniel filter is between
% 1 and 20
if ns <1 | ns>20, 
	error('number of Daniell passes (ns) must be between 1 and 20')
end



%********************* GET THE DANIELL FILTER ***********************


[spans,bwout,df,gsq]=danbw(bw,ns,tslab,padlen);
% spans specifies the spans of the filter
% bwout - the bandwidth of the estimates using the filter
% df - degrees of freedom of a spectral estimate based on the filter
% gsq - sum of squares of weights of the filter

% Call again to get the filter weights
w= danwgtn(spans);

%******************* GET PROBABILITY POINTS FOR SQUARED COHERENCY *
%
% These points depend only on the sum of squares of weights of the 
% filter used to smooth the periodogram and cross periodogram
%
p = cohprob1(gsq);


%***************  COMPUTE ROW INDICES AND OTHER INFO FOR GETTING
%		THE SEGMENTS OF THE TIME SERIES
%
% Compute some miscellaneous quantities
yrsp = yrgo + mx -1;  % ending year of time series 
tey = tslab+yrgo-1:tshift:yrsp;   % ending years of slabs
nslabs = length(tey); % number of slabs
tcy = (ceil(tslab/2)+yrgo-1):tshift:yrsp;  % central years of slabs
tcy = tcy(1:nslabs); % truncate so that all tcy are for valid
	% segments (i.e., do not over-reach ends of time series)


% Compute row-index matrix into x and y
t1 = 1:tshift:mx;  % row index to beginning years of slabs
t1=t1(1:nslabs);
t2 = t1(ones(tslab,1),:); % expand t1 to matrix
j1 = (0:tslab-1)';  % col vector
j2 = j1(:,ones(nslabs,1)); % expand j1 to matrix
T = t2+j2;  % Index to rows of x for pulling out slabs


%******************** COMPUTE FREQUENCY POINTS FOR OUTPUT
%
deltaf = 1/padlen; % frequency spacing
N = padlen/2 + 1; % number of periodogram ordinates
f = linspace(0,0.5,N); 


%********************** SIZE AND ALLOCATE FOR STORING *****

c = NaN;
C= c(ones(N,1),ones(nslabs,1));
P= c(ones(N,1),ones(nslabs,1));



%********************* LOOP OVER TIME SLABS ***************

% Compute cross-spectrum and store squared coherency and phase
for n = 1:nslabs
   t = T(:,n);  % row index in x for this slab
   zx = x(t); 
	zy= y(t);
   [cc,pp]=crspec1(zx,zy,padlen,bwout,spans,taper,0);
   C(:,n)=cc;
	P(:,n)=pp;
end

C = C .* C;  % because want squared coherency

% Check elapsed time
time1=clock-t0;


%**************** MAKE FIGURES ********************************

% Figure 1 will have the color map squared coherency
figure(1)
colormap(cool)
h1=pcolor(tcy,f',C)
[xx1 yy1]=ginput(1);
yyy1=yy1-bwout;
hl1 = line([xx1 xx1],[yy1 yyy1]);
text(xx1,yy1+0.03,'BW')
set(hl1,'Linewidth',3,'color',[0 0 0])
title(['SQUARED COHERENCY, WINDOW= ',int2str(tslab),...
	' YR,  Shift = ',int2str(tshift),' YR'])
xlabel('Year')
ylabel('Frequency (1/year)')
shading interp
hold on
Ch1=contour(tcy,f',C,p,'k');
clabel(Ch1,p)
hold off
pause


% Figure 2 will have the contour plot
figure(2)
h2=contour(tcy,f',C,p);
clabel(h2,p)
[xx1 yy1]=ginput(1);
yyy1=yy1-bwout;
hl1 = line([xx1 xx1],[yy1 yyy1]);
text(xx1,yy1+0.03,'BW')
set(hl1,'Linewidth',3,'color',[1 1 1])
title(['SQUARED COHERENCY, WINDOW= ',int2str(tslab),...
	' YR,  Shift = ',int2str(tshift),' YR'])
xlabel('Year')
ylabel('Frequency (1/year)')



% pdgm2.m

% A 386-Matlab program for univariate spectral analysis based
% on the smoothed periodogram.  Modified from pdgm1.m to give a 
% 4-frame plot of spectra for four different time series.
%
% Version 1.1, 7-23-90, by David M. Meko, Laboratory of
% Tree-ring Research, University of Arizona, 85721.
% (602) 621-2191
%
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

% ############### PRELIMINARY SETUP ####################################

% 1. Store in a ".mat" file three arrays:
%    1. A - the data array, with year as the first column
%    2. opt - the options array (see opts in variable list)
%    3. yrs - a 2-column array of first and last years for
%       each col of A (except the years column)

% 2. Start Matlab

% 3. >load  "?.mat"  , whatever you have named the ".mat" file

% 4. >pdgm2         to start program



% #################   FEATURES  ##########################################
%
% * Screen plots  of spectrum and confidence
%   intervals, with optional screen plots of intermediate series and 
%   filter weights.  
% 
% * Interactive control of degree of smoothing of spectrum and null 
%   continuum.
%
% * Optional computation of pctg of variance in selected frequency
%   ranges.
%
% * Capability of plotting any screen-viewed series to Laser-jet
%   or Epsom printer, or to HP Plotter (a general feature of 
%   386-Matlab)
%
% * Convenient selective analysis of specified columns (variables)
%   of a matrix of time series.
%   


%  ############   PRE-LOADS #########################
%
% A 	A.MAT holds the time series, each a column.  The year
% 	is held in col 1
%	

%##############   LIST OF VARIABLES  #######################

% ago	screen-prompted input value of first year of the array A
% atemp	intermediate variable used in (1) computing data window
%	(2) computing proprtion of variance in specified frequency
%	range (3) computing Daniell filters for spectral 
%	estimates, (4) computing Daniell filters for null continuum
%	(5) computing periods of spectral peaks
% btemp	intermediate variable used in (1) computing data window
% chi1	97.5 pct point of chi-sq dist.   Interpolated from lookup
%	table (see lk975.mat). Used in computing error bars around 
%	spectrum for given number of degrees of freedom. 
% chi2  	2.5 pct point of chi-sq dist...
% ctemp   intermediate variable used in (1) computing data window
% dwin	(1 x xsize-r) the data window	
% ends	(1x1i) number of values tapered on each end of X
% f	r(1 x padlen) a frequency axis for plot of DFT
% ff	r(padlen/2 x 1) frequency axis values spaced at equiv
%	Fourier freqs of padded series, and not including the
% 	symmetric other half of f.
% first	computed relative index of first year of series being analyzed
%	A "1" here would corresp to first row in A
% fline	an abscissa for plotting baseline of pdgm smoothing window
% fpeaks  r(? x 1) the frequencies of the spectral peaks
% ftemp	an abscissa for plotting inset pdgm window.  
% fyr	first year of a particular series being analyzed
% f1	r(1x1) frequency (1/yr) corresponding to p1
% f2	r(1x1) frequency... to p2
% gsq	combined variance-adjustment factor for spectral estimates.
%	Same as "g-squared" in eq 12, p. 195 in Bloomfield.
% I	i(1 x ?) the subscripts (shifted -1) of the elements of the
%	smoothed spectrum which are >= the previous element
% 	and >= the following element.  Used with yd1,yd2 and yd3 to
%	identify the peaks in the spectrum.
% i1	i(1 x ?) the subscripts of ff in the specified frequency range
%		(inclusive) in the computation of variance proportion
% igo	starting index for plot of smoothed spectrum, null line
% istop 	ending index for plot of smoothed spectrum, null line
% k	counter for loop in successively smoothing spectrum.  
%	  Incremented by one for each filter applied.
% k1	i(1 x ?) length of each successive Daniell filter applied to
%	  compute spectrum
% kill	logical variable used to control exit from while loop
%	  in smoothing periodogram and spectrum, and later in smoothing
%	  null continuum.
% l	counter for loop in successive smoothing of null continuum.
%	Similar to k for spectrum.
% last	computed index for year, last year of current series.  See first
% lk975	an 'm' file that loads a table of 97.5% points df for chi sq
% lk025	like lk975, but for 2.5% point
% l1	i(1 x ?) length of each successive Daniell filter applied to 
%	compute null continuum
% lyr	actual last year of series being analyzed.
% m	number of rows from a call to size of i1;  in call to
%	size of tab3
% mm	on first pass in loop, screen-prompted input of number 
%	of weights to smooth periodogram by to form estimated spectrum.  
%	On subsequent passes, number of weights to smooth previously
%	smoothed spectrum.
% mtap1	(1x1r) decimal fraction, percentage of series to taper.
%	mtap1/2  will be tapered on each end (Bloomfield, p. 84)
% n	number of cols from a call to size of i1;
%	  to length of y; to length of ysm1 in computing the table of 
%	  information on spectral peaks; to size(tab3) for table
% 	  of summary of spectral peaks
% ncount	counter for keeping track of how many series have been treated.
%	   Controls frame for 4-plot figure.
% nnn	i(1 x 1) which series (non-year column) of array A to do next
% opts	(3x1i) options vector, with either 1 or 0 as each element:
%	element	
%	1	screen plot intermediate series
%		  (raw data; raw data, mean subtracted; tapered; padded;
%		   unsmoothed periodogram)
%	2	compute and display % variance in specified input 
%		  wavelength range
%	3	save final screen plot as meta file for hard-copy plot
% padlen	(1x1i) padded length of X. Make a factor of 2.
% pct	r(1 x 1)  ratio of sum of pdgm ordinates specified in i1 to 
%		sum over all ordinates in y.  Equals pct variance
%		in the specified frequency range
% ppeaks	r(? x 1) periods of spectral peaks.  Computed from fpeaks.
% Pyy	c(padlen x 2) raw periodogram, using notation of Ljung
%	 Only half of Pyy contains the full info.  The other half is
%	 symmetric.
% PyyB	c(padlen x 2) raw periodogram, using Bloomfield's notation. This
%	variable is commented out, and was used originally for 
%	diagnostic purposes.  PyyB  is 1/(2*pi) times Pyy
% p1	r(1x1) screen input value of upper value of range of wavelength
% 	specified for percent variance computation.  Since this input
%	value is unlikely to match one of the discrete periodogram
%	frequencies, p1 is later converted to the period of exact
%	period of the nearest discrete pdgm frequency in the inclusive
%	specified range.
% p2	r(1x1) ... of lower value....
% straight a col vector of zeros of with length = length(wsave), used
%	to form baseline of inset pdgm smoothing window on plot
% t	r(1 x padlen) a sequential time variable, used in preparing
%	plot of DFT
% 
% tab3	a table giving fr each spectral peak:  freq, period, spectral
% 	value, null-line value, and value of lower 95% confidence
%	bound for spectrum.
% shift	a number of lags to shift the smoothed smoothed spectrum along
%	the frequency axis to compensate for the shifting implicit in 
% 	the convolution.  Needed so that smoothed spectrum plots
%	correctly.  Also used similarly for null line shifting.
% ufact  factor for adjusting variance of spectral estimate for
%	tapering data.  The original variance is
%	increased by a factor of ufact (Bloomfield, p. 194)
% var	(1x3r) variance of original series and tapered series, and
%	sum of Ljung periodogram ordinates divided by the original
%	sample size minus one.  This last element should equal the
%	variance of the tapered series.
% vdf	degrees of freedom of estimated spectrum (Bloomfield, p. 196)
% w	r(1 x mm+1) Daniel filter, same as wnew,just a shorter variable.
%	Used both in spectral estimation, and null line estimation.
% wfact sum of squares of weights in Daniell filter.  Used to 
%	adjust variance of spectral estimates. See Bloomfield (p.192).
% wnew	r(1 x mm+1) Daniell filter weights computed by convolution.
%	Used in both the estimation of spectrum and the null 
%	continuum
% wnull	r(1 x ?) Running Daniell filter for null continuum.  See wsave.
% wsave	r(1 x ?) Running Daniell filter formed by convoluting 
%	previous version of wsave wigh w.  At the end of analysis, wsave
%	will be the single filter that will yield the estimated spectrum
%	if applied to the periodogram
% X	A time series to be analyzed.  Note that this is the 
%	data series, before  removing mean, tapering, padding
% xover	the left end of the inset pdgm smoothing window will be plotted
%	at frequency xover.  Now set at 0.3.
% xsize	(1x1i) number of years of data in X, before padding   
% X1	(?x1r) tapered time series X
% X2	(padlen x 1r) padded tapered X
% y	r(padlen/2 x 1) the first half of the raw periodogram values
%	from Pyy
% ycomp	list of nul line values for spectral peaks
% yconf	listing of upper confidence level for peaks in smoothed spectrum
% 	When this value exceeds the null line for the spectrum, the 
%	specral peak is significant at 95 %
% yd1	r(padlen/2-2 x 1)  the first n-2 values of the smoothed spectrum
%	ysm1.  Used along with yd2 and yd3 to find the spectral peaks.
%	Here n is the length of ysm1, which equals padlen/2.
% yd2	r(padlen/2-2 x 1) elements 3 thru n of ysm1. See yd1.
% yd3	r(padlen/2-2 x 1) elements 2 thru n-1 of ysm1.  See yd1.
% years.dat i(nsers X 2) the beginning and ending years of valid data for
%	each series, where nsers is the total number of series.
% yext	r(3*padlen/2 x 1)  the symmetrically extended half-periodogram
%	 Formed by stacking yrev before and after y.
% yfact	0.2 times the ratio of the sizes of the maximum ordinate of the
%	ysmhi to the max ordinate of wsave.  Used to scale the pdgm
%	smoothing window for inclusion of inset on plot
% ynull	r(padlen/2+? x 1) the estimated null continuum, before shifting
% yover	amount to displace inset pdgm window in y direction for
%	iclusion on plot.  Now set to 3/4 of way to max ordinamte of
%	ysmhi
% ypeaks	r(? x1) values of spectrum at its peaks as given by ppeaks.
% yrev	r(padlen/2 x 1) the vector y with indices reversed.  Needed to 
%	stack in front of y to extend it as a symmetric vector before
%	smoothing.
% yrs	i(? x 2) first and last valid year of each time series.  The
%	vector yrs is contained in years.mat, which is loaded 
%	automatically by pdgm.m.  ? corresps to number of series in 
%	master array A, not counting the year variable.
% ysm	r(padlen/2+? x 1) the estimated spectrum.  Length depends on
%	how many and which convolutions have been done
% ysmlo	r(padlen/2 x 1) lower 95 pct confidence limit around estimate
%	spectrum in ysm1
% ysmhi	r(padlen/2 x 1) upper ...
% ysm1	r(padlen/2 x 1)  a segment of ysm covering only igo thru istop
%	Saved and used for plotting at various times.
% ysm2	r(padlen/2 x 1) a segment of ynull covering only indices igo 
%	thru istop.  Used for plot of null line.#
% ysm4	r(padlen/2 x 1) a segment of the unsmoothed pdgm used for
%	plotting.  Same as y.
% z	c(padlen x 2) DFT of X2

%##############  END OF VARIABLES LIST  ###################################


var=zeros(3,1); % to hold variances of original and tapered data, and
	% of the sum of pdgm ordinates divided by n-1, where n is
	% the original (not padded) number of years of data

hold off
clc
home
ago=input('First year in the time-series array A :');
clc
nnn=input('First series to analyze: ');


%######### TREAT EACH INDIVIDUAL TIME SERIES WITHIN A WHILE LOOP #####

mnp=220;
clg;
while nnn ~=0
ncount = 0+1;
mnp=mnp+1;

% clear up some memory space locked up in previous run through loop

clear  X i1 m n k1 l1 wnew w wsave shift igo istop yext ysm1
clear  ysm k kill l wnull ynull ysm2 ufact wfact gsq vdf chi1 chi2
clear ysmlo ysmhi ysm4 yfact xover yover straight ftemp fline yd1 yd2 yd3
clear I ff fpeaks ppeaks ypeaks yconf tab3

pause
first= yrs(nnn,1) - ago + 1;  % first year relative to first row in array
last= yrs(nnn,2) - ago +1;

fyr = yrs(nnn,1);   % actual year -- used on time-series plots
lyr = yrs(nnn,2);

disp(['Valid length of series is  ',int2str(lyr-fyr+1),' years']);
padlen=input('Desired padded length:  ');

X=A(first:last,nnn+1);


%plot(fyr:1:lyr,X)
% title('Raw Data, before tapering or padding')
 %xlabel('Observation Number')
 %ylabel('Index') 
 %text(.65,.9,['Series Number ',int2str(nnn)],'sc');
%grid
pause 
clc


% Compute variance of original series
var(1)=std(X)*std(X);

size(X);
xsize=ans(1);

mtap1 = .10;  % total decimal fraction of series to taper

% subtract mean
X=dtrend(X);
if opts(1)==1
	plot(fyr:1:lyr,X)
	title('Raw Data, Mean Subtracted');
	grid
	xlabel('Observation Number');
	ylabel('Transformed Index');
	pause
end


% Taper mtap1/2 of data on each end.  See Bloomfield, p. 84.
% Compute data window
	ends= fix((mtap1/2)*xsize);
	dwin=ones(1:xsize);
	t=1:ends;
	btemp=0.5*(1-cos((pi*(t-0.5))/ends));  %intermediate variable
	t=xsize-ends+1:xsize;
	ctemp=0.5*(1-cos((pi*(xsize-t+0.5))/ends));  % intermediate

	dwin(1:ends)=btemp;
	dwin(xsize-ends+1:xsize)=ctemp;

	clear atemp btemp ctemp

if opts(1)~=0

	 plot(dwin)
	 title('Data Window')
	 grid
	 pause 
end

% Taper by applying data window to data
	X1=X.*dwin';

if opts(1)~=0
	 plot(first:1:last,X1)
	 title('Tapered Data')
	grid 
	pause 
end


% Compute variance of tapered series
	var(2)=std(X1)*std(X1);


% pad the series to length padlen

atemp=padlen-xsize;
btemp=zeros(atemp,1);
X2=[X1;btemp];
clear X1

if opts(1)~=0
	 plot(fyr:1:fyr+padlen-1,X2)
	 title('Tapered, Padded Data')
	xlabel('Observation Number');
	 pause
end

% Compute periodogram by squaring DFT of tapered, padded series.

z=fft(X2,padlen);

t=0:1:padlen;
f=t/padlen;

if opts(1)==1
	 plot(f(1:padlen/2),z(1:padlen/2))
	 title('Discrete Fourier Tansform, Ljung Notation')
	 xlabel('Frequency')
	pause
end

% Compute the periodogram as defined by Ljung by squaring the DFT
% Bloomfield's squared dft is 1/N times Ljung's DFT.
% Bloomfield's periodogram is 1/2*pi times Ljung's periodogram.
Pyy = z.*conj(z)/padlen;


if opts(1)==1
 	plot(f(1:padlen/2),Pyy(1:padlen/2),'+g');
 	title('Unsmoothed Periodogram')
	xlabel('Frequency')
	pause
end


% compute the unsmoothed Bloomfield periodogram as 1/2pi * Pyy
% PyyB= (1.0/(2.0*pi)) .*Pyy;
% plot(f(1:padlen/2),PyyB(1:padlen/2));
% title('unsmoothed Bloomfield periodogram')


% variance of tapered (before padding) series should equal sum
% of Ljung periodogram ordinates divided by xsize, where xsize is
% the original number of years of data.

clc

if opts(2)==1  % only if you want some variance computations.
	var(3) = (1/(xsize-1)) * sum(Pyy(1:padlen));
	disp('Variance computed from data and from pdgm ordinates')
	disp('  data           pdgm')
	fprintf('\n%10g%10g',var(2),var(3))
	pause
	clc
end

% Compute fraction of variance in given range of periods

	ff=f(1:padlen/2);
	y= Pyy(1:padlen/2);  % put pdgm ordinates ...
	clear Pyy f

if opts(2)==1
	p1=input('key in longest period for variance computation:  ');
	p2=input('key in shortest period for variance computation:  ');
	clc

	f1=1.0/p1;    % frequency corresp to longest wavelength of range
	f2=1.0/p2;    % frequency corresp to shortest ...

	atemp=ones(padlen/2,1);

	i1=find(ff>=f1&ff<=f2);  %  subscripts of ff in desired freq range
	[m,n]=size(i1);          % n is number of elements found


	pct=sum(y(i1))/sum(y);  
	disp('Percentage of variance in specified period range')
	disp('      range                percentage')

	p1=1.0/ff(i1(1));
	p2=1.0/ff(i1(n));


	fprintf('\n%12f%12f%12.2f',p1,p2,pct*100)
	pause
	clc
end



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




%######  ESTIMATE THE (SMOOTHED) SPECTRUM   ###########################

% Generate a sequence of smoothed periodograms by successive application
% of Daniell filters of specified span.


k=0;
kill=0;   % will exit following while loop when kill=1

while (kill==0),
	mm=input('Span of Daniell Filter? ')	
	k=k+1;  % keep count of number of filters applied
	k1(k)=mm;   % store length of each successive filter
	
	%compute the new Daniell filter

		atemp=[.5 .5]';
	        wnew= conv(atemp,ones(mm-1,1)/(mm-1));
		% plot(w)
		% title('daniell filter')
		%pause

	% convolute Previous filter with new one, if not first in loop

		if  k==1,
			w=wnew;
			igo = n+1+(mm-1)/2;
			istop=2*n+(mm-1)/2;
			wsave=w;
		else
			w=wnew;
			shift=(mm-1)/2;
			igo=igo+shift ;
			istop=istop+shift ;
			wsave=conv(wsave,w);
		end
	
	if (k==1),
		ysm=conv(yext,w);
		
	else
		 ysm=conv(ysm,w);
	end

	%plot(ff,y,'+g',ff,ysm(igo:istop),'-r');
	%title('Smoothed Periodogram');
	%xlabel('Frequency');
	%text(.6,.9,['# of filters = ',int2str(length(k1))],'sc');
	pause
	

	if k==1
		hold on
	else
	end

	kill=input('key 1 to kill loop, 0 otherwise : ');
end
hold off
ysm1=ysm(igo:istop);
clear ysm


%##########   NULL CONTINUUM FOR SPECTRAL ESTIMATES  #######################################################

% Generate a null continuum from repeated application of Daniell
% filters until satisfied plot

l=0;
kill=0;   % will exit following while loop when kill=1
while (kill==0),
	mm=input('Span of Daniell Filter for Null Continuum ?  ')	;
	l=l+1;  % keep count of number of filters applied
	l1(l)=mm;   % store length of each successive filter
	
	%compute a Daniell filter

		atemp=[.5 .5]';
	        wnew= conv(atemp,ones(mm-1,1)/(mm-1));
		% plot(w);
		% title('daniell filter');
		%pause

% If first application of filter, save as null smoother.  If late
% application, convolve with previous filter and save as current
% filter
		
		w=wnew;
		if  l==1,
			igo = n+1+(mm-1)/2;
			istop=2*n+(mm-1)/2;
			wnull=w;
		else
			shift=(mm-1)/2;
			igo=igo+shift ;
			istop=istop+shift ;
			wnull=conv(wnull,wnew);
		end
	
	if (l==1),
		ynull=conv(yext,w);
		
	else
		ynull=conv(ynull,w);
	end

	%plot(ff,y,'+g',ff,ynull(igo:istop));
	%title('Null continuum');
	%xlabel('Frequency');
	%text(.6,.9,['# of filters so far = ',int2str(length(l1))],'sc');
	pause

	if l==1
		hold on
	else
	end


	kill=input('key 1 to kill loop, 0 otherwise : ');
end
ysm2=ynull(igo:istop);
hold off



%########## CONDIDENCE INTERVAL AROUND SPECTRAL ESTIMATES  #########################################################

% A chi-sq test is used.  Degrees of freedom depends on (1) original 
% series length,  (2) padded length, (3) proportion of data tapered,
% and (4) filter weights.

% Compute the factor gsq, that depends on the above four items.  First
% handle tapering.  See Bloomfield, p. 194.

ufact=0.5*(128-93*mtap1)/((8-5*mtap1)*(8-5*mtap1));
pause

% Variance depends on filter weights.
% Compute sum of squares of weights.

wfact=wsave'*wsave;

% Bring adjustment for padding (Bloomfield, p. 192) into computation of
% complete variance adjustment for spectral estimates (Bloomfield, p.
% 195, eq. 12).

gsq = ufact*wfact*padlen/xsize;

% Construct confidence intervals around the spectral estimates  using
% the chi-sq distribution (Bloomfield p. 196)

vdf = 2/gsq;    % equivalent degrees of freedom

%disp('Computed degrees of freedom');
%disp(vdf);

lk975;  % load table chi975 with df for 97.5 % point of chi-sq dist
lk025;  % ...        chi025....

chi1 = table1(chi975,vdf);
chi2 = table1(chi025,vdf);


ysmlo = vdf/chi1 * ysm1;
ysmhi = vdf/chi2 * ysm1;
ysm4=y;


%#########  DESIGN INSET PDGM SMOOTHING WINDOW FOR PLOT   ######################################################


% Scale the pdgm smoothing filter so can superpose on plot of
% spectrum

yfact=max(ysmhi)/max(wsave)*0.2;  % scale so that max ordinate is 1/5 
%		the size of largest ordinate of upper confidence interval

xover = 0.3;     % shift to x-ordinate of 0.3
yover= 3*max(ysmhi)/4;    % shift 3/4 of way to top for winfow plot

x1over=0.26;
y1over=1*max(ysmhi)/2;  %shift 1/2 of way to top for text plot

% Make a straight line double the length of wsave

straight=zeros(length(wsave),1);
ftemp=ff(1:length(wsave));
fline=ff(1:length(wsave));


%################  PLOT THE PDGM, SPECTRUM, CONF LIMITS, AND WINDOW ###

% next plot is without raw pdg
subplot(mnp);
plot(ff,ysm1(:),'-g',ff,ysmlo,'-r',ff,ysmhi,'-r',...
  ff,ysm2,'--w',ftemp+xover,wsave.*yfact+yover,'-w',...
  fline+xover,straight+yover,'-w')
% next is with green plusses as raw dpgm
%plot(ff,ysm4,'+g',ff,ysm1(:),'-g',ff,ysmlo,'-r',ff,ysmhi,'-r',...
 % ff,ysm2,'-w',ftemp+xover,wsave.*yfact+yover,'-w',...
 % fline+xover,straight+yover,'-w')
% title('Smoothed Periodogram with 95% Confidence Interval');
% text(.65,.85,['# Daniell filters = ',int2str(length(k1))],'sc');
text(x1over,y1over,['Region # ',int2str(nnn)]);
xlabel('Frequency (1/year)');
pause


%###########   OPTIONALLY PLOT THE FINAL DANIELL FILTER ######

	if opts(1)==1
		plot(wsave)
		title('Current Daniell Filter Weights');
		xlabel('Lag');
		pause
	end



%#############  FIND AND SUMMARIZE SPECTRAL PEAKS ################

n=length(ysm1);
yd1=ysm1(1:n-2);
yd2=ysm1(3:n);
yd3=ysm1(2:n-1);

I = find((yd3>=yd2)&(yd3>=yd1));

fpeaks=(ff(I+1))';

atemp=ones(length(fpeaks),1);
ppeaks=atemp./fpeaks;
ypeaks=ysm1(I+1);
yconf = ysmlo(I+1);
ycomp = ysm2(I+1);
tab3=[fpeaks ppeaks ypeaks ycomp,yconf];
clc
home
disp('Summary of first 17 spectral peaks');
disp('frequency, period, spectral estimate, null-line value, lower 95%')
[m,n]=size(tab3);
if m>=17, m=17;, end
disp(tab3(1:m,:));
nnn=input('next series to do--or 0 to end loop--');
end

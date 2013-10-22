% spep.m  superposed epoch analysis, with mmonte carlo confidence bands

%  Simulations based on previously modeled ARMA series

%******************   PRE-LOADS **************************************
% 
% A ... data array, first column the year
% L ... maximum lag for epoch tests
% w ... which column of A to be treated
% y ... 3x1,  first year of A, first and last years of analysis period
% k ... cv, key years
% opt ... rv, options
%	1 - npad ... # years to pad on each end of time series
%	2 - L    ... max lag for epoch analysis
%   3 - nsims ... number of simulated series to be used
%				(must be > 99)
%
%	4 - debug mode (=1 yes,  =0 no)

%*************** OTHER VARIABLES **************************
% a ... mean of x
% e ... array, random norml deviates, [length(x1) X nsims]
% g1 ... cv, real subscripts of critical prob pnts of ranked series of
%        length nsims
% g2 ... cv, diff between real subs of crit probs and next lower integ 
% gL ... cv, integer subscps just below g1
% gU ... cv, integer subspts just above g1
% k1 ... subscripts of key years relative to start of padded series
% k2 ... k1 duped to have 2*L+1 cols
% k3 ... susc array, lags added in, relative to begng of padded series
% L1 ... rv, lags -L to + L
% L2 ... array, L duped to have length(K) rows
% nlags ... 2*L+1
% p1 ... median of time series x
% P2 ... p1 converted to a rv
% TB1 ... table of means at various lags, and 5 confidence thresholds
% TB1HD ... heading for TBL1
% x ... time series to be analyzed
% x1 ... padded time series
% x2 ... cv, values of time series at various lags from each key year
% x3 ... array, x2 reshaped to length(K) x nlags
% x3M ... rv, mean value of time series relative to lag from key years
% x01,y01; x05,y05; x95,y95; x99,y99 ... stairs vars for plots ov
%   condidence bars at various prob levels

% xb, yb ... bar-plot variables for means of actual data at diff lags
% xc, yc ... stairs vars for plot of medians
% Y ... simulated time series, before adding back obs mean
% Y1 ... simulated time series, with obs mean added back

% Z1 ... cv, current simulated series within a for loop
% Z2 ... cv, values of that simulated series at each key 
   %   year and lags from the key years.  These values strung out in a 
   % single column 
% Z3 ... array, Z2 reshaped so that a row corresponds to values
   % at lags -L:L for a given key year.  Each col corresps to a lag.
% ZM ... array, mean values of each simulation at each lag. 
%		Row 1 corresponds to ave values at all lags for simulation1,
%        row 2 for simulation 2, etc.  Cols corresp to lags.
% ZM1 ... values of ZM at int subsc just below critical prob points
% ZM2 ... values of ZM at  integ subs just above crt prob points
% ZM3 ...ZM2-ZM1, diff in mean values between the two integer subscts
% ZM4 ... cutoff points of means ZM at the five critical prob levs
% ZDIFF ... incremental fraction
% ZZ ... a temporary array used for reshaping


%**********************************************************************

[m,n]=size(A);
igo= y(2)-y(1)+1;
istop=y(3)-y(1)+1;

yr=A(igo:istop,1);  % cv, years in analysis period
x = A(igo:istop,w);  % cv, time series of data, analysis pd


%********** PLOT TIME SERIES, MARKING KEY YEARS **********

a = mean (x);
plot (yr,x,'-r',yr,a(ones(length(yr),1),:));
hold on;

plot (k, a(ones(length(k),1),:),'*g');  % place marks at key years
title('TIME SERIES AND KEY YEARS');
xlabel('YEAR');
ylabel('SNOWFALL (INCHES)');
pause;
hold off;

%******** PAD SERIES AND COMPUTE YEAR INDICES ********************
%
% Pad ends of series with npad mean values.  Need to do this because
% higher lags in superposed epoch analysis will possibly extend beyond
% beginning or end of original time series.

yrbeg = y(2);  % beginning year of analysis period
npad = opt(1);
x1 = [a(ones(npad,1),:) ; x ; a(ones(npad,1),:)];

% Compute year subscripts for key years vs adjusted for padding

k1 = k' + npad(ones(length(k),1),:) - yrbeg(ones(length(k),1),:) + ...
    ones(length(k),1);

%************ FORM SUBSCRIPT ARRAY FOR KEY YEARS *******************

L = opt(2);  % max lag for epoch tests
nlags = 2* L +1
L1= [-L:L];  % rv of lags
L2 = L1 (ones(length(k),1),:);  % dupe rows of L1

k2 = k1(:,ones(nlags,1));  % dupe cols of k1

k3 = k2 + L2;  % subscr array of key years at various lags, relative to
			 % beginning of padded time series.  


%******** COMPUTE AVG AT VARIOUS LAGS FOR OBSERVED DATA ****************

x2 = x1(k3);  % data values at all subscripts given by k3.  Strung
			% out one column below the other into a single col vector.
x3 = zeros (length(k),nlags);  % to be used in reshaping
x3(:) = x2;  % reshape x2

x3M = mean (x3);  % rv, mean value of time series at each lag from key years


%*********** GENERATE SIMULATED SERIES, AND ADD BACK IN ORIG MEAN **********

% Simulated series will be of same length as padded observed series.
% ARMA model should have been estimated from approp segment of x.


nsims =opt(3);
Y=zeros(length(x1),nsims);

rand('normal');    % use normally dist random vars
e = rand(length(x1),nsims); % random normal series

for i =1:nsims
	Y(:,i)=idsim(e(:,i),th0);   % simulated time series according to specified arma
end

Y1 = Y + a(ones(length(x1),1),ones(nsims,1));  % add back obs mean


%********** AVERAGE SIMULATIONS OVER LAGS ***********************

ZM = zeros(nsims,nlags);  % pre-allocate
ZZ = zeros(length(k),nlags);  % pre-allocate

for i = 1:nsims;   % loop for each simulation
	Z1 = Y1(:,i);   % cv, current simulation
	Z2 = Z1(k3);  % values of current simulation at key years
         % and lags from key years.  
	Z3 = ZZ;  % for reshaping within for loop
	Z3 (:) = Z2;
	ZM(i,:) = mean(Z3); % means, each row corresp to a simulation, each
		% column to a lag
end

% clear ZZ Z1 Z2 Z3 Y e Y1


%**** COMPUTE CUTOFF VALUES FOR MEDIAN, 95%, 99% CONF LIMS ***********

% Aim to get a  5 by nlags array in which rows 1,2,3,4,5 corresp to
% 99 pct upper, 95 pct upper, median, 95 pct lower, 99 pct lower
% exceedance probs.  Columns correspond to lags -L thru L.  This array
% will be named ZM4.

% First rank each column of ZM.  Do this by (1) sorting ZM, 
% (2) converting the cutoff array to an array of subscripts, and
% (3) converting array of subscripts to array of cutoff values.

g1=[.01 .05 .50 .95 .99]';
g1 = g1 .* nsims(ones(5,1),:);  % cv, subscripts of critical prob
   %   points of nsim ranked values

gL = floor(g1);  % integer subscript directly below g1
gU = ceil (g1);  % integer subscript directly above g1


% ******** handle special case of less than 101 simulations

if gL(1)<1 
	gL(1) = 0;
else
end

if gU(5)>nsims 
	gU(5) = nsims
else 
end

%******* another special case

g2 = g1 - gL; % diff between real subscript and integer below

if g2(1)<0 
	g2(1)=0;
else
end

%*******  pull out relevant critical values from ZM

ZM=sort (ZM);
ZM1 = ZM(gL,:);  % mean values for int subsc just below critical prob points
ZM2 = ZM(gU,:);   % mean values for integ subs just above crt prob points
ZM3=ZM2-ZM1;   % diff in mean values between the two integer subscts
ZDIFF = g2(:,ones(nlags,1)) .* ZM3; % incremental fraction
ZM4 = ZM1 + ZDIFF;  % cutoff points for the five critical prob levs

clear g1 g2 gL gU ZM1 ZM2 ZM3 ZDIFF 


%********** CREATE BAR-GRAPH VARIABLES *****************

% Recall that x3M holds means at various lags for observed time series
%             ZM4 holds means at crit probs from monte carlo simulation
% Use median of observed as horiz base line, rather than mean

[xb,yb] = bar(L1,x3M);
p1=median(x);    % median value of observed time series
P2 = p1(:,ones(nlags,1));    % dupe p1 into a rv
yb(1)=p1;           % make takeoff points of zero the median rather than
   % the default value of zero
yb(length(yb))=p1;
[xc,yc] = stairs (L1,P2);
  yc(1)=p1;
  yc(length(yc))=p1;

[x01,y01] = stairs(L1,ZM4(1,:));
[x05,y05] = stairs(L1,ZM4(2,:));
[x95,y95] = stairs(L1,ZM4(4,:));
[x99,y99] = stairs(L1,ZM4(5,:));

y01(1,[1 length(y01)])=p1(:,ones(1,2));
y05(1,[1 length(y05)])=p1(:,ones(1,2));
y95(1,[1 length(y95)])=p1(:,ones(1,2));
y99(1,[1 length(y99)])=p1(:,ones(1,2));

%*********  MAKE TABLE OF MEANS AND CONF LEVELS ***********

TB1HD = ['   LAG       MEAN      .01      .05      .50     .95      .99'];
TB1 = [L1' x3M' ZM4'];


%********* MENU FOR RESULTS *******************************

j2=1;
while j2~=4
	j2 = menu('Select an Item','Bar Plot with 95% CL',...
         'Bar Plot with 99% CL',...
		'Table With Sample Means and Confidence Bands',...
		'All Done');
	
	if j2==1
		plot(xb,yb,'-r',x05,y05,'--g',x95,y95,'--g',xc,yc,'-r');
		title('MEANS AT VARIOUS LAGS FROM KEY YEARS');
		xlabel('LAG');
		ylabel('MEAN GROWTH INDEX');
		pause;
		clg
	elseif j2==2
		
		plot(xb,yb,'-r',x01,y01,'--g',x99,y99,'--g',xc,yc,'-r');
		title('MEANS AT VARIOUS LAGS FROM KEY YEARS');
		xlabel('LAG');
		ylabel('MEAN GROWTH INDEX');
		pause;
		clg
	elseif j2==3
		clc
		home
		disp(TB1HD);
		disp(' ');
		disp(TB1);
		pause
	

	elseif j2==4
		disp('You Done')
	
	end
end

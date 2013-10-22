function resin(z,th,M,maxsize)
%RESID  Computes and tests the residuals associated with a model
%
%	E = resid(Z,TH)
%
%	Z : The output-input data Z=[y u], with y and u being column vectors.
%	For multi-variable systems Z=[y1 y2 .. yp u1 u2 ... un]. 
%	For time-series Z=y only.
%	TH: The model to be evaluated on the given data set. (Format as
%	described by HELP THETA)
%	E : The residuals associated with TH and Z. [resid(Z,TH); just performs
%	and displays the tests, without returning any data.]
%
%	The autocorrelation function of E and the cross correlation between
%	E and the input(s) is computed and displayed. 99 % confidence limits 
%	for these values are also given (based on the hypothesis
%	that the residuals are white and independent of the inputs). These 
%	functions are given up to lag 25, which can be changed to M by
%	E = resid(Z,TH,M). The correlation information can be saved and re-
%	plotted by
%	[E,R] = resid(Z,TH).         resid(R);
%
%	E = resid(Z,TH,M,MAXSIZE) changes the memory variable MAXSIZE from
%	its default value. See HELP AUXVAR.

%	L. Ljung 10-1-86,1-25-92
%	Copyright (c) 1986-92 by the MathWorks, Inc.
%	All Rights Reserved.

% ** Set up the default values **

[E,R1]=resid(z,th);

[NN,nz]=size(z);
N1=R1(1,26);

if (NN~=N1),
	error('N from r1,r2,r3 must match sample length in z')
end


% *** Compute the residuals and the covariance functions ***

[mrr1,nrr1]=size(R1);
M=nrr1-2;  % will plot lags 0 thru M

% ** Confidence interval for the autocorrelation function **
%
% R1 as returned by resid.m holds covariances, so must scale by 
% dividing by variance
% 99% limit on acf is 2.58 divided by sqrt(sample length)

nr=0:M-1; % x axis for plotting
cl99=  2.58 / sqrt(N1);
cl99 = cl99(ones(M,1),:);  % 99% cl as col vector


r1=R1(:,1:M);


% Get acf of residual out of r1.  r1(1,1) holds the
% variance of residuals, while r1(1,:) holds the autocovariances

a1 = r1(1,:)/ r1(1,1);  % acf of resids


% Plot acf of residuals the 99% confidence bands
subplot(2,1,1)
stem(nr,a1)
hold on
plot(nr,[cl99 -cl99],'y:')
%axes('Position',[0.1 0.1  0.65  0.4])
%plot(nr,a1','--',nr,a2','-',nr,a3',':',nr,[cl99 -cl99],':r')
%legend('OE(1,0)','OE(1,1)','OE(2,1)');
title ('RESIDUALS ANALYSIS -- CONF LIMITS 99%')
xlabel('Lag (Yr)');
ylabel('Autocorrelation')
hold off


%*********   cross-correlation plot **********************

nr=-M+1:M-1;  % x-axis for ccf plots

% Row index into r1 pointing to acv of noise, 
% ccv (ind2 and ind3) and acv of input
ind1=3;  ind2=2;  indy= 1;  indu=4;
 
% get the relevant ccfs out of r1
c1= [r1(3,M:-1:1) r1(2,2:M)]/(sqrt(r1(indy,1) * r1(indu,1)));


% *** Compute confidence lines for the cross-covariance functions **


sdreu=2.58*sqrt(r1(indy,1)*r1(indu,1)+2*(r1(indy,2:M)*r1(indu,2:M)'))/sqrt(N1)*ones(2*M-1,1);
cc199= sdreu / (sqrt(r1(indy,1)*r1(indu,1)));

subplot(2,1,2)
stem(nr,c1')
hold on
plot(nr,[cc199 -cc199],'y:');
ylabel('Crosscorrelation with Input')
xlabel('Lag(Yr)')
hold off


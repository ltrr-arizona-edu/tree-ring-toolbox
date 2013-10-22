% smpd.m  computes and smoothes a periodogram
% requires input from screen:
%   Assumes tsa of data in in X
%   Assumes filter is in A


X=[2 5 3 1 7 3  3  1  3  6];
A=[.2 .2 .2  .2 .2];
SX=size(X)
SA=size(A)


% Filter shorter than time series?

if SX(2)>SA(2)
	B=conv(X,A)
else
	disp('Hey Jack, filter is longer than time series.')
	pause
end


% Filter has odd number of weights

if rem(SA(2),2)==0
	disp('Hey Jack, filter must have odd number of weights')
	pause 
else
end


% Sum of filter weights must equal 1.0

asum=sum(A)

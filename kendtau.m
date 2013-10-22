function [T,tau] = kendtau(x,y)

% Kendall's tau.  See Practical Nonparametric Statistics, second
% edition, by W. J. Conover, p. 256.   
%
%  D. Meko 12-12-93
%
%
%****************   INPUT ARGS ************
%
% x, y -- column vectors, each a time series covering same interval
%
%
%***************************  OUTPUT ARGS ********
%
% T (1 x 1)   number of concordant minus number of discordant pairs. T is
%	the statistic entered in table A12 of Conover to test for significane
% tau (1 x 1) kendall' tau
%
%


n=length(x);
if length(y) ~= n,
	error('x and y must be same length col vectors');
end

[r,j]= sort (x);  % sort one variable from smallest to largest.
%   j holds corresponding indexes into original x

s= y(j);  % elements of y rearranged into order of sorted x


% Form logical matrix that will allow pointing out 'ties' in ranked x
R1= r(:,ones(n,1));  % dupe cv r into matrix
R2=R1';
RI= R1==R2;

% Form logical matrices that will allow pointing out of concordant and
% discordant pairs
s1=s'; % s to row vector
S1= s1(ones(n,1),:);   % dupe rv s1 into matrix
S2=S1';

% Identify concordant and discordant pairs, minus 'ties'
SC = S1 < S2 & ~ RI;  % concordant
SD = S1 > S2 &  ~ RI;  % discordant

% Convert upper right triangular parts of SC, SD to zeros so that
% concordant or discordant pairs are summed over succeeding observations
SC=tril(SC,-1);
SD=tril(SD,-1);

% Sum the concordant pairs; sum the discordant pairs
nc = sum(sum(SC));
nd = sum(sum(SD));

% Compute the "T" value that is tested in table A12 of Conover
T = (nc-nd);

% Compute the tau statistic
tau = T / (n* (n-1)/2);


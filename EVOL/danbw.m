function [spans,bw2,df,sos]=danbw(bw1,ns,nobs,padlen)
% Spans and properties of a Daniell filter that has ns equal spans, and
% a sum-of-squares of filter weights close to that of a rectangular
% filter with a bandwidth of bw
%
% D Meko 7-10-95
%
%************** IN ARGS *********************
%
% bw1 (1 x 1)r   ideal desired bandwidth, in "per-year" frequency units

% ns (1 x 1)i   number of passes of the Daniell filter
% nobs (1 x 1)i original (unpadded) length of the time series
% padlen (1 x 1)i	length of the series under consideration; typically
%				the length after padding
%
%******************  OUT ARGS *********************
%
% scans (1 x ns) lengths of the individual Daniell filters to
%	be convoluted. Note that all elements of scans will be equal,
%	and are constrained to be odd by this function. Thus [3 3 3 3]
%	might be a returned value of scans
% bw2 (1 x 1)r bandwidth of selected Daniel filter
% df  (1 x 1)r degrees of freedom of Daniel estimate
% sos  (1 x 1)r sum of squares of Daniel filter weights
%

%
%************  NOTES ***************************************
%
% If a series has length n, the Fourier frequencies are spaced at 1/n,
% and a rectangular filter with bandwidth bw covers a known number n1
% of these frequency points. This rectangular filter has a sum of 
% squares of weights of 1/n1.  The ideal Daniel filter is defined as 
% one giving the same (or nearly the same) sum of squares
%
% Consider the following variables:
%
% ns - number of passes of the Daniell filter (specified)
% n1 - length of the rectangular filter, which is a minimum for
%		possible nonzero span of the desired final Daniel filter
% n2 - length of the individual ns spans of the Daniel filter (unknown)
% nmin  - length of Daniell filter than convoluted ns times has same
%    nonzero span as the n1-length rectangular filter

% Compute number of (padded) Fourier frequencies covered by the desired bw for
% a rectangular filter
n1 = round(bw1/(1/padlen)); 

% Compute the length of the Daniell filter that convoluted ns times
% would give the same non-zero span

nstart = (n1+ns-1)/ns;
nstart=floor(nstart);

% Make sure length of filter is odd
if rem(nstart,2)==0;
	if nstart>1
		nstart=nstart-1;
	else
		nstart = 1;
	end
end


sos1 = 1/n1; % equals sum of squares of n1-length rect filter
sos2 = sos1+1;; % starting value for sum of squares of daniell filter
sos2old = sos2;
n2 = nstart-2;

k=1;
while (k); % while sos of daniell filter greater than rect
	n2=n2+2;
	sos2old = sos2;
	f=danwgtn(n2(:,ones(ns,1)));
	sos2 = sum(f .* f);
	if sos2<sos1
		k=0;
	end
end

diff1 = abs(sos2-sos1);
diff2 = abs(sos2old-sos1);
sos = sos2;

if diff2 <diff1;
	n2 = n2-2;
	sos=sos2old;
end

% Bandwidth of final Daniel filter
n5=1/sos; % Length of hypoth rect filter with same sum of squares 
	% as selected Daniel filter
bw2 = n5/padlen;  % Equivalent bandwidth of the Daniel filter

% Degrees of freedom of the Daniel filter
df = (n5*2)* (nobs/padlen);



spans = n2(:,ones(ns,1));




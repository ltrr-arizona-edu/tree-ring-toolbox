function f = danwgtn(spans)
%
% filter weights for a modified-Daniell filter, single or multi-span
%
% D Meko 7-7-95
%
%******************** IN ARGS ********************************
%
%	spans (1 x ns)i spans of the individual Daniell filters that
%		will be convoluted to get final filter.  Each element of 
%		spans is the number of weights in a single-span filter.
%		Elements of spans must be odd and positive.
%
%****************** OUT ARGS *********************************
%
% f(n x 1)r the filter weights of the resultant Daniell
%		filter
%
%************ USER-WRITTEN FUNCTIONS CALLED
%
% danwgt1.m -- weights of single-span daniell filter




[dum,ns]=size(spans);
if dum~=1, error('spans must be a row vector'), end;
if any(spans<1)  |  any(rem(spans,2)==0)
	error('Elements of spans must be positive and odd');
end

% Handle trivial case of unit filter
if all(spans==1)
	f=1.0;
	return
end

% Compute half-lengths of individual filters
mvect = (spans-1)/2;


% Compute  initial filter
if spans(1)==1; % First span is equal to 1
	f=1.0;
else
	f = danwgt1(spans(1));
end

% Loop over remaining elements of spans, if necessary
if ns ==1;
	return
else
	for k = 1:(ns-1);
		f1=f;
		s2 = spans(k+1); % 
		if s2==1
			f2=1.0;
		else
			f2=danwgt1(spans(k+1)); 
		end
		f = conv(f1,f2);
	end
end



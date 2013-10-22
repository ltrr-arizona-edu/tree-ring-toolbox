function f = danwgt1(n)
%
% filter weights for a single-span modified-Daniell filter
%
% D Meko 7-7-95
%
%******************** IN ARGS ********************************
%
%	n (1 x 1)i span (length) of filter. Must be odd and positive.
%		Note that the half-length is (n-1)/2,
%
%****************** OUT ARGS *********************************
%
% f (n x 1)r the filter weights
%
%

% Check args
if n<=0 | rem(n,2)==0,  
	error('Span must be odd and positive'), 
end


% Initialize
f = zeros(n,1); 

wt = 1 / (n-1); % weights except for the two end weights
f=[ wt/2; wt(ones(n-2,1),:); wt/2];






function y = moddan(x,sym,spans)
%
% Filter a series with a modified-Daniel filter
%
% D Meko 7-7-95
%
%************************** IN ARGS ***********************
%
% x (mx x 1)r  series before filtering
% sym (1 x 1)i symmetry for series extension before filtering
% 		1=even, 0=zeros, -1=odd [see extend.m]
% spans (1 x ns)i  spans of Daniel filters to be convoluted to build
%		final filter.  Individual spans must be positive and odd.
%
%****************  OUT ARGS **********************************
%
% y (my x 1)r  smoothed series.  Same length as x.
%
%
%********************* USER-WRITTEN FUNCTIONS CALLED *****
%
% danwghtn.m  -- daniell filter
% extend.m   -- extend ends of series before filtering


% Check and size
[mx,nx]=size(x);
if nx ~=1 | mx <=1,
	error('x must be a col vector')
end
[ms,ns]=size(spans);
if ms~=1, error('spans must be a row vector'), end;


% Trivial case of all elements of spans unity -- unit filter
if all(spans==1);
	y=x;
	return
end

% Get the Daniell filter weights
f = danwgtn(spans);
n =length(f); % number of weights in filter
m = (n-1)/2;  % half-length of filter

% Extend x before filtering it
x1 = extend(x,sym,m);

% Convolute filter with extended series
x2 = conv(f,x1);
nx2 = length(x2);

% Cull out valid part of x2 
y = x2(n:nx2-n+1);





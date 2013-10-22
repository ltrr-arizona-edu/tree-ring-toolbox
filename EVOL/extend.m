function y=extend(x,sym,mext);
% Extend a series (say, a fft, or a time series) 
%
% D Meko 7-7-95
%
%****************  IN ARGS *************
%
%	x (mx x 1)r  series to be extended
%	sym (1 x 1)i type of symmetry.  1=positive, -1=negative,
%		0=extend with zeros
%	mext (1 x 1)i number of values to extend each end of x by
%
%
%********************** OUT ARGS ****************
%
%	y (my x 1)r  extended series.  Equivalent to x with mext values
%		tacked on each end
%
%**************  NOTES ************************************
%
% Even symmetry: end values reflected across x ordinate of last value .
% Odd symmetry: end values reflected across y ordinate of last value,
% then across xordinate of last value.
%
% Tip for use of extend.m with odd-number-of-weights Daniell filter
% Say the time series to be filtered is x, the filter has nf weights,
% and the output smoothed series is to w.  Do this:
%
%  1. y = extend(x,1,(nf-1)/2) % extend series by half-length of filter
%  2. z = conv(f,y);  % convolute filter with extended series
%  3. w = z(nf:(mz-nf+1)); % cull out valid values of z



a=NaN;
[mx,nx]=size(x);
if nx~=1, error('x must be col vector'), end
if mx <= mext, error('length of x must exceed mext'), end
if sym <-1 | sym>1
	error('sym must be -1,0, or 1')
end

% Set a few row indices
m1 = mx + 2*mext; % length of extended series
y = a(ones(m1,1),:);  % initialize extended series
m1head = 1:mext; % row index in y of leading extended values
m1tail = (m1 -mext+1):m1; % row index of trailing extended values 
m1body = (mext+1):(m1-mext);  % row index of unchanged values


% Build extended series
if sym ==0; % extend with zeroes
	y = [zeros(mext,1); x; zeros(mext,1)];
elseif sym ==1; % even symmetry
	y(m1head)=flipud(x(2:(mext+1)));
	y(m1tail)=flipud(x((mx-mext):(mx-1)));
	y(m1body)=x;
elseif sym==-1; % odd symmetry
	goseg = flipud(x(2:(mext+1))-x(1));
	endseg = flipud(x((mx-mext):(mx-1)) - x(mx));
	y(m1head)= x(1) - goseg;
	y(m1tail) = x(mx) - endseg;
	y(m1body)=x;
end

if any(isnan(y));
	error('A NaN has been found in extended series y')
end


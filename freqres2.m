function R=freqres2(w,deltat,f)
% freqres2: frequency response of a symmetrical discrete-valued smoothing or filtering function
% R=freqres2(w,deltat,f);
% Last revised 10-12-00
%
% Frequency response of a symmetrical discrete-valued smoothing or filterin function.
%
%*** INPUT
%
% w (1 x mw)r filter weights.  Number (mw) of weights must be odd. Weights must be
%    symmetrical and discete valued. (e.g.:  0.25 0.50 0.25)
% deltat (1 x 1)r  the data interval, or time between successive observations (e.g. 1 year)
% f (nf x 1)r frequencies (cycles/time unit) at which frequency response is to be computed
%
%
%*** OUTPUT
%
% R (nf x 1)r amplitude of frequency response at frequencies f
%
%
%*** REFERENCES 
%
% Panofsky, H. and Brier, G., 1958.  Some applications of statistics to meterorology.
% The Pennsylvania State University.

% Check that number of weights odd
[nw,mw]=size(w);
if nw~=1;
   error('w must be row vector');
end;
if rem((mw+1),2) ~= 0;
   error('Number of weights must be odd');
end;

% Check that weights sum to 1
if abs(sum(w)-1) > 1e-10;
   error('weights must sum to 1');
end;

  


% Check frequencies
nf=length(f);
if any(f<0) | any(f>0.5);
   error('f must be between 0 and 0.5');
end;


% Check that weights symmetrical
ic = (mw+1)/2; % index of central weight
iw = (ic+1):mw; % index of weights to right of center
iother = fliplr(1:(ic-1));
wcheck = w(iother);
wc = w(ic);
w = w(iw); 
d = abs(w-wcheck);
L=(d<=1e-8);
if ~all(L)
   error('w not symmetrical');
end;

n=length(w); % number of weights on one side, not counting center

% Matrix of index to weights
k=(1:n)';
K = repmat(k,1,nf);

% Matrix of frequency
F = repmat(f',n,1);

% Matrix of weights
W = repmat(w',1,nf);

% Part to be summed
B = W .* cos((K .* F)*2*pi*deltat);
[mB,nB]=size(B);

if mB==1;
   R=wc + 2*B;
else;
   R = wc + 2*sum(B); % a rv of freq response
end;

R= R';







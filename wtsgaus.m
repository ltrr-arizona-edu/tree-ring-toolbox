function b=wtsgaus(p)
% wtsgaus: weights for gaussian filter with 50% response at period of p years
% b=wtsgaus(p);
% Last revised 10-13-00
%
% Coefficients, or weights, of binomial filter of length n
%
%*** INPUT
%
% p (1 x 1)i  period (years) at which filter is to have amp frequency response of 0.5
%
%*** OUTPUT
%
% b (1 x n)r  computed weights
% 
%
%*** REFERENCES
% 
% WMO 1966, p. 47
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- stats
%
%*** NOTES
%
% Amplitude of frequency response drops to 0.50 at a wavelength of 
% about 6 standard deviations of the appropriate guassian curve


% Check that period of 50% response at least 3 yr
if p<3;
   error('Period of 50% response must be at least 3 yr');
end;

sigma=p/6;  % Gaussian curve should have this standard deviation

x=-100:100;
b=normpdf(x/sigma,0,1);
bmax=max(b);
bkeep = b>=0.05*bmax; % keep weights at least 5% as big as central weight
b=b(bkeep);
b=b/sum(b); % force weights to sum to one


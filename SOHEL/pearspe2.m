function [r,npsp] = pearspe2(x,y,k)
%
% USAGE : [r,npsp] = pearspe2(x,y,k)
%   Pearson product-moment and Spearman correlation coefficient
%   between paired columns of a time series matrix.
%
%
% INPUTS 
%-------
%       x( mx x nx )    - a time series matrix
%       y( my x ny )    - a time series matrix, same size as x
%       k( 1  x 1  )    - Option
%                           1 = Pearson
%                           2 = Spearman
%
%
% OUTPUTS
%--------
%       r( nx x 1 )     - correlation coefficients
%                           r(1) is between column 1 of x and y,
%                           r(2) is between column 2 of x and y, etc.
%
%	npsp(1 x 4 )	- A vector indicating 4 significance levels
%			  npsp(:,1) -> 20 % S.L. in 2-tailed test
%			  npsp(:,2) -> 10 % S.L. in 2-tailed test
%			  npsp(:,3) -> 5 % S.L. in 2-tailed test
%			  npsp(:,4) -> 1 % S.L. in 2-tailed test
%
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% PEARSP.TAB	Table containing the critical values
%____________________________________________________________________

  if( nargin ~= 3 );
    error('Wrong number of input arguments.');
  end

  if( size(x) ~= size(y) );
    error('Size of matrix x and y are not the same.');
  end
  
  [mx,nx] = size(x);
  if( k==2 );                   % Must convert x,y to ranks.
    
    [s,x] = sort(x);
    [s,y] = sort(y);

    mat = (1:mx)';
    mat = mat(:,ones(1,nx));
    matstring(:) = mat;

    a = (0:nx-1);
    a = a(ones(mx,1),:);
    a = a * mx;

    x = x + a;
    t1(:) = x;
    t11(t1) = matstring;
    x = zeros(mx,nx);
    x(:) = t11;

    y = y+a;
    t1(:) = y;
    t11(t1) = matstring;
    y = zeros(mx,nx);
    y(:) = t11;

  end

  xmean = mean(x);
  ymean = mean(y);
  xmean = xmean(ones(mx,1),:);
  ymean = ymean(ones(mx,1),:);
    
  r = sum( (x - xmean) .* (y - ymean) ) ./ (std(x) .* std(y));
  r = r / (mx-1);
  r = r';
  
end                            

if ~exist('pearsp'),
  load pearsp.tab;
end

npsp=table1(pearsp,mx);

% End of the function

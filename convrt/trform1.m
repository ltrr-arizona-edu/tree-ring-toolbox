function Y=trform1(X,khow,a,b);
% trform1:  log10, inverse log-10, or power transformation of a variable
% Y=trform1(X,khow,a,b);
% Last revised 9-23-00
%
% Log-10, inverse log-10, or power transform of a time series
%
%*** INPUT
%
% X (mX x nX)r  matrix of time series to be transformed;  can be a col vector
% khow(1 x 1)i  option for transformation
%    ==1 y = log10(x+a); % log10 transform
%    ==2 y = exp(log(10)* x))-a; % inverse log10 transform
%    ==2 y =  (x+a)*exp(b);  % power of b transform
% a (1 x 1)r  <==0> number to be added to x before log10 transform or raising to power
% b (1 x 1)r  (<==1>; moot if khow=1): power that (x-a) is to be raised to 
%
%*** OUTPUT
%
% Y (mY x nY)r  matrix of transformed time series
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES

% Check X dimensions
[mX,nX]=size(X);
if mX==1; 
   error('X must be matrix or col vector');
end;
nser = nX; % number of series
nobs = mX; % number of observations

if all([khow~=1 khow~=2 khow~=3]);
   error('know must be be 1, 2 or 3');
end;


if khow==1;  % LOG 10 TRANSFORM
   b=1;
   X1=X+a; % add increment to each series (usually a==0);
   if any(X1<=0);
      error('X1 must be >0 for log10 transform');
   end;
   Y = log10(X);
elseif khow==2;% inverse LOG 10 TRANSFORM
   b=1;
   Y = exp(log(10)*X)-a;
elseif khow==3;
   % If use negative power in transform, data multiplied by -1 after
   % raising to power to avoid reversing direction of relationships
   Xa = X+a;
   % To avoid odd behavior in power-transformed series, do not allow the
   % oringal series (with a added) to have both positive and negative values.
   % Can imagine with power of 2 transform that -2 and +2 both transformed
   % to +4
   % I simplify this by requiring that all (X+a)>=0.  
   if any(any(Xa<0));
      error('All elements of (X-a) must be >0 for power transform');
   end;
   if b==0; % raising to power of zero give constant time series at 1
      error('Not allowed to raise to power of zero');
   elseif b>0;
      Y = (X+a) .^ b;
   else;
      if any(X+a)==0;  % zero raised to neg power is infinite
         error('zero (X + a) yield infinite value for b negative');
      else;
         Y = -(X+a) .^b; % negate to avoid reversing direction of anomalies
      end;
      
      
   end;
   
end;



   
   
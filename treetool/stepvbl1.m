function inext = stepvbl1(y,X,Lin,Lmask)
% stepvbl1: choose next variable to enter regression based on maximum R-squared 
% inext = stepvbl1(y,X,Lin,Lmask);
% Last revised: 3-15-98
%
% Low-level fuction called by other functions (e.g., respfun1.m).  One possible
% use is in a forward-stepwise regression calling function to decide which
% predictor variable to enter next.
% 
%
%*** IN **************
%
% y (my x 1)r    predictand
% X (mX x nX)r   tsm of potential predictors;  note that my must equal mX
% Lin (1 x nX)L   predictors already in model (Lin=1, otherwise Lin=0)
% Lmask (1 x nX)L  potential predictors to be masked out of consideration as
%   next variable to enter equation
%
%*** OUT ******************
%
% inext (1 x 1)i  index to column of X pointing to the variable chosen as
%   next to enter the regression
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- none
%*** TOOBOXES NEEDED -- statistics
%
%*** NOTES ************
%
% The number of columns in X is the total number of variables.  Lin has
% that same number of rows, with 1 meaning the variable is in the model and
% 0 that the variable is not in.
%
% Variables marked as 0 in Lin are candidates for entry as the next variable, except
% tha Lmask can indicate not to consider a variable for entry. Lmask entry of 1
% means do not allow the corresponding X variable to enter the equation
%
% stepvbl1 calls MATLAB function  regress.m to compute the regression R-squared


if nargin~=4;
   error('nargin should return 4');
end

%-----------  Input checks

% X and y
[mX,nX]=size(X);
[my,ny]=size(y);
if ny~=1; 
   error('y must be a col vector')
end
if mX ~= my;
   error('X and y must be same row size');
end

% Lin and Lmask
if ~islogical(Lin) | ~islogical(Lmask);
   error('Lin and Lmask must be logical');
end
if size(Lin,2)~=nX | size(Lmask,2)~=nX;
   error('row size of Lin and Lmask must equal col size of X');
end
if size(Lin,1)~=1 | size(Lmask,1)~=1;
   error('Lin and Lmask must be row vectors');
end


% Check that no NaN in X or y
if any(isnan(y));
   error('y contains NaNs');
end
if any(any(isnan(X)));
   error('X contains NaNs');
end


% Check that have at least one variable available after allowing for 
% those in the equation (Lin==1) and those masked out (Lmask==1)
Lavail = ~Lin & ~Lmask;  
if ~any(Lavail);
   error('No variables available after allowing for those in and those masked');
else
   navail = sum(Lavail);
end

% If only one variable available, that's the one to pick
if navail==1;
   inext = find(Lavail);
   return;
end

% Size vector to save R-squared values of candidate models
s  = repmat(NaN,navail,1);  

if ~any(Lin);
   X1=[];
else
   X1 = X(:,Lin);  % Store the subset of variables already in model
end


X1=[ones(my,1) X1];  % put a ones vector on the predictor matrix

iavail = find(Lavail);  % 

% Loop over the potential new variables
for n = 1:navail;
   itry = iavail(n);
   xtry = X(:,itry); % potential new variable
   X2=[X1 xtry];
   
   [b,rint,r,rint,STATS] = regress(y,X2,0.01);
   
   s(n) = STATS(1);  % Store R-squared
end

[ybig, ibig]=max(s); % maximum value of the R-squared

inext = iavail(ibig);
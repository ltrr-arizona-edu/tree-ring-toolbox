function yn = rescaly(x,y)
%
% USAGE : yn = rescaly(x,y)
%   Rescales y to have the same mean and standard deviation 
%   as x
%
%
% INPUTS
%-------
% x (nx x 1)	Vector of SKE Data
% y (ny x 1)	Vector of IND or RW data (ny = nx)
%
%
% OUTPUTS
%--------
% yn (ny x 1)	Rescaled output 
%
%
% NO USER WRITTEN FUNCTIONS NEEDED
%_______________________________________________ 

xmn = mean(x);
ymn = mean(y);

yn = (y - ymn)/std(y);
yn = (yn * std(x)) + xmn;

% End of file
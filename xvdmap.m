% xvdmap.m  computes an array with following columns, where each
%  row corresps to one of the 238 xvalid points:

%  1. x coordinate
%  2. y coordinate
%  3. mean estimation error (mean for 275 years)
%  4. mean squared estimation error
%  5. standard error of the estimate 
%  6. std deviation of observed (actual) data
%  7. ratio of "5" to "6"

%****************   PRELOADS *****************************

% xvd.mat - data array
% xvdxy.dat - x, y coords of the points, in cols 1,2

%**********************************************************

X=zeros(238,7);    % preallocate


X(:,1:2)=xvdxy (:,1:2);
X(:,3)=ksu(:,4) /275.0;
X(:,4)=ksq(:,4) /275.0;
X(:,5)=sqrt(X(:,4));

a=ksq(:,1)/275;  % mean square of observed data at each point
b=(ksu(:,1)/275) .^2; % squared mean of observed data  

X(:,6)=a-b;  % variance equals mean square minus squared mean

X(:,6) = sqrt(X(:,6));  % standard deviation of observed data

X(:,7)= X(:,5) ./ X(:,6); % ratio of std error of estimate to 
 % standard deviation of observed

plot(1:238,X(:,3));
title('MEAN ESTIMATION ERROR AT 238 POINT');
pause;

plot(1:238,X(:,4));
title('MEAN SQUARED ESTIMATION ERROR AT 238 POINTS');
pause


plot(1:238,X(:,5));
title('STANDARD ERROR OF THE ESTIMATE AT 238 POINTS');
pause

plot(1:238,X(:,6));
title('STANDARD DEVIATION OF OBSERVED DATA AT 238 POINTS');
pause

plot (1:238,X(:,7));
title('RATIO, STD ERROR EST TO STD DEV OBSERVED');

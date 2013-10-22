function [j,r,n]=highestr(R,ord,i)
% highestr.m:  subfunction of surrgt2.m.  Identify surrugate time series from correlations
% CALL:  [j,r,n]=highestr(R,ord,i);
%
%********************* IN *****************************
%
% R (? x 4)r  mtx of correlation coef, sample size, series1 pointer, series2 pointer
%   as returned from r4sov.m
%
% ord (? x 1)i;  "order" of time series pointed to by cols 3,4 of R.
% i (1 x 1)i:  pointer to the key series
%
%
%********************** OUT ****************************
%
% j (1 x 1)i  pointer indicating the highest correlated surrogate series
% r (1 x 1)r  the correlation coefficient between key series and the selected surrogate
% n (1 x 1)i  sample size (yr) for r
%
%***************** NOTES ******************************************
%
% Surrogate is the series of next lowest order to that of series i with the
% highest correlation with series i

[mR,nR]=size(R);

% Compute number of time series
nsers = length(ord);

% Get the order of key series
keyorder=ord(i);
if keyorder == 0;
   error('keyorder must be greater than zero');
elseif keyorder>3;
   error('keyorder must be less than 4');
end

% Get the order of series for comparison
surorder = keyorder-1;

% Any series of lower order?
i1 = ord<=surorder;
if sum(i1)==0,
   error('No series of lower order than key series');
end

nsurr = sum(i1); % number of potential surrogate series

% If only one series of lower order, bingo, that is the surrogate
if nsurr==1;
   j=find(i1);
   % get pointer to correct row of R
   L1= (R(:,3)==i | R(:,4)==i) & (R(:,3)==j | R(:,4)==j) & ~isnan(R(:,1));
   r=R(L1,1);
   n=R(L1,2);
   return
else; % more than one lower order series
   
   % Find which correlations (rows of R) apply to the key and potential surrogate series
   j=find(i1); % row vector of pointers to potential surrogate series
   J=j(ones(mR,1),:); % expand the rv to a matrix, same row size as R
   % Make separate matrices, same size as J, of cols 3 and 4 of R
   r3=R(:,3);
   R3=r3(:,ones(nsurr,1));
   r4=R(:,4);
   R4=r4(:,ones(nsurr,1));
   % Logical pointer to rows of R applying to any of the surrogate series
   L1=(any((J==R4)'))'  |   (any((J==R3)'))';
   % Logical pointer to rows of R applying to the key series
   I=i(ones(mR,1),:); % dupe pointer to key series into cv, same row size as R
   L2= r4==i | r3==i;
   % Logical pointer to rows of R with correlations between key and a surrogate
   L3=L1  & L2;
   % Get cv subset of correlation coefficients for key vs potential surrogates
   rs = R(L3,1);
   % Get corresponding row subset of R
   Rsub=R(L3,:);
   
   % Find maximum R
   [r,imax]=max(rs); % sorts in ascending order; rmax is the correlation,
   % imax is its position in the original vector rs
      
   % Get the pointer to the selected surrogate series
   rsel=Rsub(imax,3:4); % 1 x 2 rv; a pointer to key and the other series
   j = rsel(rsel~=i);  % pointer to the selected surrogate
   n= Rsub(imax,2); % sample size (# yr) for r
      
   
   
   
end


function I = rowtsm01(YRS)
% rowtsm01: compute row indices of series in a stacked tsm
% CALL: I = rowtsm01(YRS);
%
% Meko 2-23-98
%
%*************  IN 
%
% YRS (? x 2)i or (? x 3)i   start and end years of each series in the stacked matrix
%     A third row indicating start row index sometimes exists
%
%*********** OUT 
%
% I (? x 2)i start and end row index of each time series 
%
%
%*********** NOTES 
%
% User has a matrix with time series stacked one after another. These might be
% monthly precip records for several stations, for example.  YRS has the start
% and end year for each series.  It is assumed that no missing years (rows) 
% are in the referenced data matrix of stacked time series.
%
% User wants row index matrix that will allow him to easily pull out the data
% for each series

%----------  Compute number of series
[m1,n1]=size(YRS);
nser = m1;

%------------ YRS can have 2 or 3 cols
if  ~(n1==2 | n1==3);
   error('YRS must have 2 or 3 cols');
end

% Check that ending year no earlier than start year
X = YRS(:,1:2)';
if any(diff(X)<0);
   error('YRS cols 1,2 inconsistent as start and end years');
end




% Computations
nyr = YRS(:,2)-YRS(:,1)+1;% number of rows for each time series
endi = cumsum(nyr);  % row index of last year for each series
begi = endi - nyr +1;  % row index of first year for each series

I = [begi endi];

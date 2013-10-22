function [R,N]=simptch(X,nmin,khow)
% simptch:  matrix of similarity statistics between pairs of time seris with "patchy" coverage
% R=simptch(X);
% Last revised 4-25-01
%
% Given a time series matrix of core indices.  The indices cover various periods, so that the matrix has NaNs in
% places.  The overlap period for pairs of series varies, and some pairs may not overlap at all.  You want the statistics 
% summarizing the similarity of variations in pairs of series.  Some minimun length of overlap is required for
% a statistic to be computed and stored.  Alternative statistics: correlation coefficient, median cross-product of
% departures from 1.0, relative departure (homegrown), coefficient of departure (homegrown). 
%
%*** INPUT
%
% X (mX x nX)r   time series matrix, NaN filled where no data
% nmin (1 x 1)i  minimum allowable sample size for a valid statistic of similarity between series
% khow (1 x 1)i  option for statistic
%   ==1 correlation coefficient
%   ==2 median cross-product of departures from 1.0 (see notes)
%   ==3 relative departure (see notes)
%   ==4 coefficient of departure (see notes)
%
%
%*** OUTPUT
%
% R (nX x nX)r  correlation coefficients, or NaN.  R(i,j) holds correlation between series i and j
% N (nX x nX)i  sample size for the correls in R;  NaN if no valide
%
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%
%*** NOTES
%
% Statistics, assuming series are x and y 
% Correlation coefficient, R
%   R is the Pearson product-moment correlation, computed in the usual way.  
% Cross-product of departures, CP 
%   CP is the sum of the products of (x-1.0) and (y-1.0) for all overlap points
% Relative departure, RD 
%   Similar to the "reduction of error" statistic. Compute mean square departure of series x from series y as MSD1.
%   Compute mean square departure of the indiviual values of x and y from 1.0 as MSD2.  Compute RD as 
%   1-(MSD1/MSD2)
% Coefficient of departure, CD
%   Similar to coefficient of correlation, except that the numerator is the sum of products of 
%   departures from 1.0 rather than from the sample means.  As with the correlation coefficient, the
%   denominator scales the statistic by the product of the standard deviations.
%
%
% Aim.  Function written to deal with similarity in smoothed variations of tree-ring indices.  Normal for 
% tree-ring indices is 1.0.  So two smoothed series that might have slightly opposing trends but both be way below 
% 1.0 are actually "similar", despite posssibly having negative product-moment correlation coefficient.  The 
% alternative statistics covered by khow are an attempt to gauge the similarity.



% Size & check input
[mX,nX]=size(X);
[mtemp,ntemp]=size(khow);
if ~(mtemp==1 & ntemp==1);
    error('khow must be scalar');
else;
    if khow<1 | khow>4;
        error('khow must be between 1 and 4');
    end;
end;
if nmin<10;
    error('nmin must exceed 10');
end;

% Allocate
N = repmat(NaN,nX,nX);
R= repmat(NaN,nX,nX);



% LOOP OVER SERIES

for n = 1:nX; % Loop over series
    x = X(:,n); % current key series
    L1 = ~isnan(x); % mark valid data
    nsum1=sum(L1); % full non-NaN length of series x
    if nsum1<nmin; % if fewer than min allowable sample size for this series
        % no action needed; R (m,:) remains NaN
    else;
        x=x(L1);
        % Pull Slab of tsm corresponding to non-NaN rows of X, for cols beginning next series above key series
        Y = X(L1,((n+1):nX));
        [mY,nY]=size(Y);
        % Punt if no paired series as adequate sample overlap
        L2 = ~isnan(Y);
        if   all(sum(L2)<nmin); 
            % No action needed; R(m,:) remains NaN
        else; % At least one paired series has enough overlap for computation of a valid correlation
            for m = 1:nY; % loop over series higher in sequence number than key series
                y = Y(:,m);
                L2this = L2(:,m);
                nsum2 = sum(L2this);
                if nsum2<nmin; % sample size of overlap too small
                    % No action needed, R(n,m) remains NaN
                else; % Can get a valid correlations
                    nsum = nsum2; % sample size for the correlation
                    xthis = x(L2this); % the overlapping part of the key series
                    ythis = y(L2this); % the overlapping part of other series
                    % Departures from 1
                    xd1 = xthis-1.0;
                    yd1 = ythis-1.0;
                    
                    if khow==1; % correlation coeff
                        r = corrcoef(xthis,ythis);
                        r = r(1,2);  % the correlation coefficient between series
                        
                    elseif khow==2; % median cross product of departures from 1.0
                        % Scaling  by 100 makes scale of result closer to range of -1 to 1 
                        % (e.g., indices of 1.1 and 1.1 give 100*(.1 * .1) = 1 rather than 0.01
                        r =    100*median(xd1 .* yd1);
                    elseif khow==3; % relative difference
                        MSD1 = mean((xthis - ythis)  .^2); % mean squared departure of series x from series y
                        D = [xd1; yd1]; % concatenate departures from 1.0 for x and y series
                        MSD2 = mean(D .^ 2); % mean squared departure of all data values from 1.0
                        r = 1 - (MSD1/MSD2);  % mean square departure statistic
                    elseif khow==4; % coefficient of departure from 1.0
                        cp = mean(yd1 .* xd1); % mean product of departures from 1.0
                        r = cp / (std(xthis) * std(ythis)); % scaled by standard devs
                    end;
                    R(n,(n+m))=r;
                    N(n,(n+m))=nsum;
                    
                    
                end;
            end; % for m=1:nY
        end; % all(sum(L2)<nmin);
    end; % nsum1<nmin
end; % n = 1:nX
    
                    
                





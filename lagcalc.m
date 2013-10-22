function kj=lagcalc(kopt,nvars,nneg,npos,ijk)
% lagcalc: given a matrix of lagyr3.m-made lagged predictors, compute lag from column
% CALL: kj=lagcalc(kopt,nvars,nneg,npos,ijk);
%
% Meko 10-7-97
%
%***************** IN ***********************************
%
% kopt(1 x 1)i option for mode
%   kopt(1)==1; compute col number in lagged mtx given lag and variable number
%          ==2; compute lag given variable number and col in lagged matrix
% nvars (1 x 1)i number of variables (sites) in the lagged matrix.  This is 
%   number of zero-lagged variables
% nneg (1 x1)i  number of negative lags in lagged matrix (e.g., 2);
% npos (1 x 1)i  number of positive lags in lagged matrix (e.g., 2);
% ijk (1 x 3)i  ijk(1) -- variable number (possible range 1 to nvars)
%               ijk(2) -- col number in lagged matrix 
%               ijk(3) -- lag
%
%*************** OUT *****************************
%
% kj (1 x 1)i  output whose definition depends on kopt(1)
%   kopt(1)==1: kj is the computed col number in lagged predictor matrix
%   kopt(1)==2: kj is the lag of the variable whose column in the lagged
%         predictor matrix is given
%
%************** NOTES *********************************8
%
% lagcalc.m written as auxillary function for use with reglag1.m in cross-
% referencing columns in a lagged predictor matrix to lags on particular tree-ring
% sites
%
% Example of use in kopt(1)==1 mode:
%  kj = lagcalc(1,57,2,2,[3 NaN -1])--> Lagged predictor matrix was made up of
%     lags -2,-1,0,1,2 on 57 tree-ring sites, for total col size of 57*5.  You
%     want to know the column in that lagged matrix holding the lag t-1 data for
%     site 3.  Call would return kj==60, meaning column 60
%
% Example of use in kopt(1)==2 mode:
% kj = lagcalc(2,57,2,2,[3 60 NaN]) --> Same setup of lagged predictor matrix.
%     You have the column number in the predictor matrix as col 60, and know
%     that your site is site 3.  What is the lag represented by col 60?
%     Call would return -1, meaning lag t-1 years
%
% lagcalc.m assumes that lagged predictor matrix has been built as in lagyr3.m, which
% means for n sites, lag-0 variables are in cols 1-n,  lag -1 variables in the next
% n cols, lag -2 variables in the next n cols, etc.... then the lag +1, lag +2, etc
% variables

% 
%
%
%
[m1,n1]=size(ijk);
if m1~=1 | n1~=3;
   error('ijk must be 1 x 3');
end



i=ijk(1);  % site number
j=ijk(2);  % col number (is NaN if kopt==1)
k=ijk(3);  % lag (may be negative, zero,  positive, or NaN (if kopt==2)

if i<=0 | i>nvars;
   error('Illegal site number -- inconsistent with nvars');
end
if j<=0 | j>nvars*(1 +nneg+npos);
   error('Illegal col number -- inconsistent with nvars');
end
if k<0; % if neg lag
   if abs(k)>nneg;
      error('k inconsistent with nneg');
   end
end
if k>0;
   if k>npos;
      error('positive k greater than npos');
   end
end


%********************************************8


if kopt(1)==1; % mode 1:  given site number and lag, compute col in lagged matrix
   if ~isnan(j);
      error('For kopt(1) equal to 1, input j should be NaN');
   end
   if k==0; % zero-lag year
      j=i;  % column
   elseif k<0; % negative lag
      j = i + (-k) * nvars; % column
   else; % positive lag
      j = i + (nneg*nvars) + k*nvars;
   end
   kj=j;
elseif kopt(1)==2;  % mode 2: given site number and col number in lagged matrix, compute lag
   if ~isnan(k);
      error('For kopt(1) equal to 2, input k should be NaN');
   end
   if j<=nvars; % must be lag zero
      if i~=j;
         error('input i must equal input j if j less than nvars');
      end
      k = 0;  % lag is zero
   elseif j>nvars & j<=nvars * (nneg+1);  % lag is negative
      ktemp = -(j-i)/nvars;
      if rem(ktemp,1)~=0;
         error('Specified col j inconsistent with site number i');
      else
         k=ktemp;
      end
      if abs(k)>nneg;
         error('Computed k less than negative nneg');
      end
      
   elseif j>(nvars * (nneg+1));  % lag is positive
      ktemp = (j - i - (nneg*nvars))/nvars;
      if rem(ktemp,1)~=0;
         error('Specified col j inconsistent with site number i');
      else
         k=ktemp;
      end
      if k>npos;
         error('Computed k greater than npos');
      end
      
   end
   kj=k;
end

   

      
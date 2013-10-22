function [Ix,Iy]=indseq1(mx,my,igox,igoy,ncal,minoff)
% indseq1:  index-sequential row-index matrices for bivariate data analysis
% CALL: [Ix,Iy]=indseq1(mx,my,igox,igoy,ncal,minoff);
% Last revised: 12-26-99
%
% A utility function used to build row indices to time series.  The row indices can 
% be used by the calling function to pull index-sequential samples from the series, 
% for example, to estimate error bars around coefficients for response functions.
%
%*** IN
%
% mx (1 x 1)i  number of observations in first time series in calling function -- x
% my (1 x 1)i  number of obs in second ...  -- y
% igox (1 x 1)i  row of first year of calibration period in x
% igoy (1 x 1)i row of first year of calibartion period in y
% ncal (1 x 1)i  number of years in calibration period
% minoff (1 x 1)i  minimum offset of any index-sequential sample from the calibration 
%   period observed juxtaposition of x,y.  Typically minoff==1
%
%*** OUT
%
% Ix (ncal x nsamp)i  each col of Ix is row-index to an index-sequential sample of x
% Iy (ncal x nsamp)i  each col of Iy is row-index to an index-sequential sample of y
%   Matching columns of Ix,Iy for a single sample for drawing an index-sequential pair
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Need.  Wanted to compute confidence limits for response function analysis (respfun1.m)
% by the circular-index-sequential method.  Say have a climate series y and a tree-ring 
% series x.  Say a response function is calibrated for some ncal-year subset of overlap
% of x and y.  The index sequential error bars are computed (1) making x and y circular 
% by linking the last year to the first year, and (2) repeatedly shifting x relative to 
% y and doing the calibration for all possible non-observed juxtapositions.  
%
% minoff.  This argument included to handle possible highly autocorrelated series such 
% that index-sequential samples at small offsets from observed are not too similar to 
% the observed.  In practice, best to simply set minoff==1 and accept consequences if 
% if a slight shift in x relative to y gives just as good a relationship as no shift.
%
% Algorithm makes sure no element of any x,y index-sequential sample pair has the same 
% offset as the observed data.  This ensures that index-sequential samples are not 
% contaminated by "true" matchup of x and y.


%*** INPUT CHECK *****************

if length([mx my igox igoy ncal minoff]) ~=6;
   error('All input args should be scalar integers');
end;

if mx<5 | my<5;
   error('Minimum acceptable mx or my is 5');
end;

if ncal<5;
   error('Minimum acceptable ncal is 5');
end

if mx-igox+1 <ncal;
   error('ncal must be at least mx-igox+1');
end;

if my-igoy+1 <ncal;
   error('ncal must be at least my-igoy+1');
end;

if minoff<1 | minoff>5;
   error('minoff cannot be greater than 5 or less than 1');
end;


% Build some needed column vectors
x=(1:mx)';  y=(1:my)';   % index vectors same size as x and y
x = [x;x]; y=[y; y];  % stack x and y to double length
ix = (1:(2*mx))'; iy=(1:(2*my))';  % row indices to double-length x, y


% Build col vector of possible starting row-indices for index-seq samples of x
i1 = 1:mx;
I1=repmat(i1,my,1);
i1=I1(:);

% Likewise for y
i2=(1:my)';
i2=repmat(i2,mx,1);

% Use the starting indices i1, i2 and calibration length ncal to build row indices 
% to index-sequential samples of length ncal
idope = (0:(ncal-1))';
Iinc  = repmat(idope,1,(mx*my));
I1=repmat(i1',ncal,1)+Iinc;
I2=repmat(i2',ncal,1)+Iinc;
% ... but note that elements of I1 may range up to 2*mx and of I2 to 2*my.

% Make references in I1,I2 fold back
L1 = I1>mx;
if any(any(L1));
   I1(L1)=I1(L1)-mx;
end;
L2 = I2>my;
if any(any(L2));
   I2(L2)=I2(L2)-my;
end;
Ix=I1;
Iy=I2;

% Delete as unacceptable index sequential samples any that contain an x and y 
% element with the same offset as the observed data
offset = igox-igoy;
ix = Ix(1,:); iy = Iy(1,:);
A=ix-iy  - offset;
B=(ix-mx)-iy  - offset;
C=ix - (iy-my)   - offset;
LA=abs(A)<minoff;
LB=abs(B)<minoff;
LC=abs(C)<minoff;
L=LA | LB | LC;
if any(L);
   Ix(:,L)=[];
   Iy(:,L)=[];
end;

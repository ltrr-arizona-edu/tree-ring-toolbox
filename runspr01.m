function [y,yry]=runspr01(x,s,yrx,m,xc,pmin,kopt) 
% runspr01:  runs probability for a reconstruction, using reconstruction error 
% [y,yry]=runspr01(x,s,yrx,m,xc,pmin,kopt);
% Last revised 11-30-99
%
% For a given run-length, computes probability that any year is the ending year
% of that many consecutive years outside the threshold
%
%*** INPUT
%
% x (mx x 1)r  reconstructed time series
% s (mx x 1)r  standard deviation of errors (e.g., RMSE of validation)
% yrx (mx x 1)i  year vector for x and s
% m (1 x 1)i  length of run of interest (e.g., 6)
% xc (1 x 1)r  critical value (threshold) for defining a single-year event
% pmin (1 x 1)r minimum probability to accept the year as possibly being part
%   of an m-year event (e.g., pmin==0.05)
% kopt (1 x 1)i  options
%   kopt(1) whether event is value below or above threshold
%     ==1 below
%     ==2 above
%
%*** OUTPUT
%
% y (my x 1)r  prob that year yry is end of m consecutive values below xc
% yry (my x 1)r  year vector for y
%
%*** REFERENCES - NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Example of use.  Have reconstructed flow series. See that worst modern drought 
% was 6 consecutive low-flow years -- below 12 MAF. Ask if such 6-year drought is 
% likely to have happened before.  
%
% Error standard deviation, s.  Function samples a normal distribution with a given 
% mean and standard deviation to find the probability that any given year was below 
% the threshold xc.  You might have various types of estimates of s. For example, 
% RMSE of cross-validation, or standard error of the estimate.
%
% m -- length of run of interest.  Output series y will have mx-1 fewer values than 
% input time serie x because mth value of x is the first that could possibly be 
% identified as the last year of an m-year event.
%
% Definining m-year event.  The probability that m consecutive years were dry is computed 
% as the product of the individual probabilities that each of the m years was dry.
%
% pmin -- You must specify some value (minimum or max, depending on kopt(1)) that you 
% require the reconstructed value to be outside before the year will even be considered 
% part of an m-year event.  If the probability for all m years in any m-year period 
% does not reach at least pmin, the m-year event probability is set to zero for that 
% period.


% QUALITY CONTROL

[mx,nx]=size(x);
[ms,ns]=size(s);
[mtemp,ntemp]=size(yrx);
if nx~=1 | nx~=1 | ntemp~=1;
   error('x, s and yrx  must be column vectors');
end;
if mx~=ms | mx ~=mtemp;
   error('x, s and yrx must be same length');
end;
if any(isnan(x)) | any(isnan(s)) | any(isnan(yrx));
   error('x, s and yrx not allowed to have NaNs');
end;

clear mtemp ntemp;
[mtemp,ntemp]=size(m);
if mtemp~=1 | ntemp~=1;
   error(' m must be scalar');
end;
[mtemp,ntemp]=size(xc);
if mtemp~=1 | ntemp~=1;
   error(' xc must be scalar');
end;
[mtemp,ntemp]=size(pmin);
if mtemp~=1 | ntemp~=1;
   error(' pmin must be scalar');
end;

if m<1 | m> mx/2; 
   error('m must be at least 1, and must be less than half the length of x');
end;
if xc>max(x) | xc<min(x);
   error('xc cannot be greater or less than all values in x');
end;
if pmin<0 | pmin>1;
   error('pmin is a probability, must be between 0 and 1');
end;

[mtemp,ntemp]=size(kopt);
if mtemp~=1 | ntemp~=1;
   error('kopt must be row vector of length 1');
end;
if ~(kopt(1)==1 | kopt(1)==2);
   error('kopt(1) must be 1 or 2');
end;


% COMPUTE SINGLE-YEAR PROBABILITIES
p1 = normcdf(xc,x,s);
% p1 is a column vector, same length as x, with the prob that the year is a drought year

% Build lagged matrix of single-year probs
P = repmat(NaN,mx-m+1,m);
for j = 1:(m);
   jgo = m-j+1;
   jsp = mx-j+1;
   P(:,j)=p1(jgo:jsp);
end;

yry=yrx;
yry(1:(m-1))=[];


% BUILD VECTOR OF PROB THAT ANY GIVEN 6-YEAR PERIOD WAS DROUGHT EVERY YEAR

% Replace any P smaller than pmin with 0
L = P<pmin;
if any(any(L));
   P(L)=0;
end;
clear L;

% Product across the m years for every ending year
y = (prod(P'))';
   

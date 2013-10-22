function [r,df,pvalue,sigmar]=corrts1(x,y,kopt,kdir)
% corrts1:  significance test of correlation between two time series
% [r,df,pvalue,sigmar]=corrts1(x,y,kopt,kdir);
% Last revised 7-2-01
%
% Test of null hypothesis that population correlation is zero against either one-sided or
% two sided alternative hypothesis.  Optional adjustment of degrees of freedom for significance
% tet for first order autocorrelation in the individual series
%
%*** INPUT
%
% x (mx x 1)r one time series
% y (my x 1)r the other time series, same length as x
% kopt(1 x 2)i
%   kopt(1) ?-tailed test
%       ==1 one tailed
%       ==2 two tailed
%   kopt(2) adjust df for autocorrelation?
%       ==1 yes
%       ==2 no
% kdir (1 x 1)i if one sided, which direction for alter hypoth
%   ==1  population r greater than 0
%   ==2  population r less than 0
%   Note: set kdir==[] for two sided test 
%
%
%*** OUTPUT
%
% r (1 x 1)r correlation coefficient
% df (1 x 1)r  degrees of freedom for test of significance of r
% pvalue(1 x 1)r  p-value for test 
% sigmar (1 x 1)r computed standard deviation of theoretical correlation coefficient
%
%*** REFERENCES
%
% Dawdy, D.R., and Matalas, N.C., 1964, Statistical and probability analysis of hydrologic data, part III: 
% Analysis of variance, covariance and time series, in Ven Te Chow, ed., Handbook of applied hydrology, 
% a compendium of water-resources technology: New York, McGraw-Hill Book Company, p. 8.68-8.90.
%
% Panofsky, H.A., and Brier, G.W., 1958, Some applications of statistics to meteorology: 
% The Pennsylvania State University Press, 224 p.
%
%*** UW FUNCTIONS CALLED -- NONE
% acf.m -- autocorrelation function
%
%*** TOOLBOXES NEEDED
% Statistics
%
%
%*** NOTES
% 
% Null hypothesis is that the samples are uncorrelated.  Test valid only for this null hypothesis. 

% CHECK INPUTS

% x
[mx,nx]=size(x);
if nx~=1 | mx<5;;
    error(' x must be cv of minimum length 5');
end;
nsize1 = mx;

% y
[my,ny]=size(y);
if nx~=1 | mx~=my;
    error(' y must be cv of same length as x');
end;
if any(isnan(x)) | any(isnan(y));
    error('NaNs in x or y');
end;


% kopt
if size(kopt,1)~=1 | size(kopt,2)~=2;
    error('kopt must be 1 x 2');
end;
if ~all(kopt>=1 & kopt<=2)'
    error ('kopt entries must be 1 or 2');
end;


% Compute Correlation
r =corrcoef([x y]);
r=r(1,2);


% Adjust Sample size for autocorrelation?
if kopt(2)==1; % adjust
    % HERE CODE CALL TO ACF ETC
else;
    N =nsize1;
end;

% Compute standard deviation of pop correl coef
df = N - 2; % degrees of freedom
sigmar = sqrt(1/df);
meanr=0; % theoretical mean of r is zero


% Tails 
if kopt(1)==1;
    ntail=1;
else;
    ntail=2;
end;

pcdf = normcdf(r,meanr,sigmar);

if ntail==2; % two sided test
    if pcdf>0.5;
        d=(1-pcdf)*2;
    else;
        d=pcdf*2;
    end;
elseif ntail ==1; % One tailed
    
    if kdir==1; % alt that r>0
        d=1-pcdf;
    else;
        d=pcdf;
    end;
    
end;

pvalue=d;









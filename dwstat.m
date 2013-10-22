function [d,choice]=dwstat(e,durbin,k)
% dwstat:  Durbin-Watson statistic on residuals of a MLR model
% CALL: [d,choice]=dwstat(e,durbin,k);
%
%*****************  IN *******************
%
% e (nobs x 1)r  regression residuals
% durbin (38 x 11)r  table of .01 (one sided) values for D-W stat
% k (1 x 1)i  number of predictors in the regression model
%
%****************** OUT *******************
%
% d (1 x 1)r  computed D-W statistic
% choice (1 x 1)s   decision on the null hypothesis that rho==0
%   'R' reject
%   'U' uncertain
%   'A' accept
%
%****************** NOTES *******
%
% The alpha level is interpted as 0.01 for a one-sided test and 0.02 for a two-sided test
%
% dwstat.m handles only up to 5 predictors.  This is because that's the limit of the
% table I keyed in from Draper and Smith
%
% The table for interpreting the D-W statistic is assumed to be loaded by the
% calling function.  For reference, this table is durbin.dat in the structure 
% durbin, stored in \mlb\tables



[nobs,ntemp]=size(e);
if ntemp~=1;
   error('e must be col vector');
end

if nobs<15 | nobs>100;
   error('My Durbin-Watson table only covers sample sizes 15 to 100');
end
if k>5;
   error('My Durbin-Watson table only covers up to 5 predictors in model');
end



%************** Compute D-W stat (see Ostrom 1978, p 27, eq 2.24)

dtop = diff(e);  % first difference of residuals
num =    sum(dtop .* dtop); % sum of squares of them
denom = sum(e .* e); % sum of squares of resids
d = num/denom; % D-W stat


%************* GET APPROPRIATE COLUMNS OF TABLE

klow = k*2;
khi = klow+1;

ylow = durbin(:,klow);
yhi = durbin(:,khi);

dlow = interp1(durbin(:,1),ylow,nobs);
dhi = interp1(durbin(:,1),yhi,nobs);

if d<dlow ;
   choice='R';
elseif d>= dlow & d<dhi;
   choice='U';
elseif d>=dhi & d<= (4-dhi);
   choice = 'A';
elseif d>(4-dhi) & d<(4-dlow);
   choice='U';
elseif d>(4-dlow);
   choice='R';
end

function [C,p]=chisqdst(x,dist,kgroup,nbin,c,cinc);
% chisqdst: chi-squared test for agreement of sample distribution with theoretical
% [C,p]=chisqdst(x,dist,kgroup,nbin,c,cinc);
% Last revised  7-6-99
%
% Tests for agreement of sampling and theoretical distributions using a chi-square test.
% You can choose from normal or lognormal distribution for the comparison.  You can 
% specify the number of bins for the test, and how those bins are set up.  Theoretical
% distribution fit using parameters of the sample (e.g., mean and std dev for normal
% distribution).
%
%*** INPUT
% 
% x (mx x 1)r  observations
% dist(1 x ?)s  theoretical distribution.  Options:
%   'normal', 'lognormal'
% kgroup (1 x 1)i  option for how bins are selected
%   ==1  specify upper limits of each bin (c).  Values
%        < c and not in a lower bin fall in that bin.  Highest bin infinity.
%        Set nbins, uinc==[].
%   ==2  specify number of bins (nbins), upper limit for lowest bin (u), and 
%        increment of upper limit to form other bins (cinc). 
%   ==3  specify number of bins (nbins).  Evenly spaced percentiles of fitted 
%        theoretical distribution are used as upper class limits. Thus, nbins==10
%        makes the limits at percentiles 10,20 ... 100.  c and cinc set to [].
%   ==4  specify number of bins (nbins).  nbins even-spaced bins are set up
%        between the minimum and maximum values of the sample.  c and cinc set to [].
% nbins (1 x 1)i  number of bins. Might be [] depending on kgroup.
% c (1 x 1)r or (nc x 1)r  upper class limit for lowest class, or for all classes,
%      depending on kgroup
% cinc (1 x 1)r  increment for upper class limits.  Might be [], depending on kgroup.
%
%*** OUTPUT
%
% C (1 x 1)r  chi square value
% p (1 x 1)r  alpha-value for C



%---- Fit theoretical distribution
if strcmt(dist,'normal');
   [uhat,sigmahat,u,u]=normfit(x);
else;
   error('normal only dist allowed');
end

%--- Make upper class limits
if kgroup==4; % base on sample min and max; make evenly spaced
   if ~isempty(c) & ~isempty(cinc);
      error('c and cinc must be empty for kgroup==4')
   end;
   xmin = min(x); xmax=max(x);
   upper = linspace(xmin,xmax+1,nbins+1);
   upper(1)=[];
elseif kgroup==1;  % specify upper limits of bins
   if ~isempty(nbins) & ~isempty(cinc);
      error('nbins and cinc must be empty for kgroup==1')
   end;
   upper = c;
   nbins=length(c);
elseif kgroup==2;  % specify upper limt of lowest bin and increment of limits
   b=[o:(nbins-1)]';
   upper = c  +  b*cinc;
elseif kgroup==3;  % quantiles of theoretical distrib fit to sample
   if strcmp(dist,'normal');
      q=linspace(0,(nbins+1),1); % quantiles for upper limits 
      upper = norminv(q,uhat,sigmahat);
      upper(1)=[];
      upper(nbins)=upper(nbins)+0.01*upper(nbins); 
   end;
end;

   
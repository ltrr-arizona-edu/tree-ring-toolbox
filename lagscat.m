function lagscat(x,m,nfigs)
% lagscat: lagged scatterplots for a single time series
% lagscat(x,m,nfigs);
% Last revised  9-14-99
%
% Lagged scatterplots are scatterplots of x(t) vs x(t-k). Plots are made 
% for lags k=1,...m.  m may set as high as N/3, where N is the sample size
%
%*** IN ********************
%
% x (m1 x 1)r  time series, m1 observations
% m (1 x 1)i  number of maximum lags for scatterplots 
% nfigs (1 x 1)i <optional>  start new figure windows at window nfigs+1
%
%*** OUT ********************
%
% No output arguments
% Plots are created in windows nfigs+1, nfigs+2,....
% Four pairs of axes are produced in each figure window
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% More on the input:
%  m: maximum allowable value is 13
%
%  nfigs:  allows you to retain nfigs existing figure windows when calling 
%  lagscat from another function that has generated graphics you do not want 
%  to overwrite.  For example, if calling function has 6 
%  figure windows you want to retain, set nfigs=6 in calling lagscat.  
%
%  Approximate 95% confidence level for correlation coefficient is computed as 
%  2/sqrt(N), where N is the sample size (e.g., number of years)



%*******************  CHECK INPUTS

% If only 2 input args, assume you want scatterplots to start in figure 1
if nargin==2; 
   nfigs=0;
end;

[m1,n1]=size(x);
if n1~= 1;
   error('x should be a column vector');
end
nsize  = m1; % number of observations

% m must be <= 1/3 the size of x
if m>m1/3;
   error('m not allowed to be greater than 1/3 the sample length');
end;



% Find the maximum and minimum of the series so that can set the 
% x and y axis limits same for each figure
xmax = max(x);
xmin = min (x);


%********************** BUILD LAGGED TIME SERIES MATRIX
kwind=0;
for n  = 1:m;
   if rem((n-1),4)==0;
      kwind=kwind+1;
      k = 1;
      figure(kwind+nfigs); % open figure window, with fig number adjusted by nfigs
   else
      k=k+1;
   end
   
   % Get series segments
   u = x(1:(nsize-n));   
   v = x((n+1):nsize);
   
   nyr = length(u) ;  % number observations for the scatterplot
   strn = sprintf('%5.0f',nyr);
   
   
   % Compute correl coef
   rpear = corrcoef([u v]);
   rpear = rpear(1,2);
   strr = sprintf('%4.2f',rpear);
   
   % Compute 95% sign r
   r95 = 2/sqrt(nyr);
   str95 = sprintf('%4.2f',r95);
   
   % Build string for title
   str1 = ['\itr = \rm' strr ',  \itN = \rm' strn ',  \itr_{95} = \rm' str95];
   
   
   subplot(2,2,k);
   plot(u,v,'.');   
   lsline; % plot best fit straight line
   txtx = ['x(t-' int2str(n) ')'];
   txty = 'x(t)';
   title(str1,'FontSize',8);
   
   xlabel(txtx);
   ylabel(txty);
   
   % Set axis limits
   set(gca,'Xlim',[xmin xmax],'Ylim',[xmin xmax]);
   
   % Make the plot box square instead of rectangular
   set(gca,'PlotBoxAspectRatio',[1 1 1]);

  
end

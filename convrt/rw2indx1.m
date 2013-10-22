function [Y,S]=rw2indx1(X,kopt,w);
% rw2indx1: matrix of ring widths to matrix of tree-ring indices
% Y=rw2indx1(X,kopt,w);
% Last revised 1-25-02 
% 
% Time series matrix of ring widths converted to matrix of indices by optional straight-line or spline
% detrending, and optional ratio or difference method.  Needed as utility function by function maxwave2.m,
% which investigates the covariation of indices at various frequency ranges
%
%*** INPUT
%
% X (mX x nX)r  time series matrix of ring width measurements
% kopt (1 x 2)i  options  for detrending
%   kopt(1) functional detrending
%       ==1 horizontal line at ring-width mean
%       ==2 least squares straight line
%       ==3 spline of set wavelength (see w)
%   kopt(2) ratio vs difference
%       ==1 ratio
%       ==2 difference (see notes)
% w (1 x 1)r  wavelength (yr) of spline;[] if kopt(1)==1
%
%
%*** OUTPUT 
%
% Y (mY x nY)r  time series matrix of indices, same size as X
% S (2 x nY)r   means (row 1) and standard devs (row 2) of the indices in Y

%*** NOTES
%
% Difference detrending.  If kopt(2==2, index is shifted and scaled too. Shifted to mean of 1.0.
%  scaled to have save variance as if detrended by ratio method using same trend line
%
% w -- the wavelength at which the frequency response of the spline has amplitude 0.5


%---- CHECK INPUT 

[mX,nX]=size(X);
if mX<20; 
    error('Row size of X must be at least 20');
end;

for n = 1:nX;
    x = X(:,n);
    % Check for internal NaN
    k=intnan(x);
    if k==1;
        error(['Internal NaN in col ' int2str(n) ' of X']);
    end;
    % Check for number of values in series
    ngood = sum(~isnan(x));
    if ngood<20;
        error(['Fewer than 20 nonNaN values in col ' int2str(n) ' of X']);
    end;
end;

if kopt(1)==3; 
    if isempty(w);
        error('w must be specified if kopt(1)==3');
    end;
    if w<6; 
        error('w <6 gives unstable cubic smoothing spline');
    end;
end;


%---- COMPUTE THE TREND LINES

G = repmat(NaN,mX,nX); % to hold trend line

if kopt(1)==1; % horizontal line (null detrending
     for n = 1:nX;
        x=X(:,n);
        L=~isnan(x);
        x=x(L);
        mx = length(x);
        xh = mean(x);
        G(L,n)=xh ; % store trend line
    end;

elseif kopt(1)==2; % straight line method
    for n = 1:nX;
        x=X(:,n);
        L=~isnan(x);
        x=x(L);
        mx = length(x);
        t=(1:mx)';
        U=[ones(mx,1) t];
        b=U\x;
        xh=b(1)+b(2)*t; % trend line
        G(L,n)=xh ; % store trend line
    end;
    
elseif kopt(1)==3; % spline method
    p=splinep(w,0.5); % spline parameter
    for n = 1:nX;
        x=X(:,n);
        L=~isnan(x);
        x=x(L);
        mx = length(x);
        t=(1:mx)';
        xh=csaps(t,x,p,t);
        G(L,n)=xh ; % store trend line
    end;
  
end;

%---- CONVERT MATRIX OF RW TO INDICES, BOTH BY RATIO AND DIFF METHODS

Hr=  X ./ G;
Hd = X - G;
[mH,nH]=size(Hr);

% Compute the means and statdard deviations of indices
Sr = nanstd(Hr);
Sd = nanstd(Hd);
Mnr = nanmean(Hr);
Mnd = nanmean(Hd);


%---- SCALE INDICES, IF NECESSARY, IF DIFFERENCE METHOD

if kopt(2)==1;
    Y=Hr;
    S=[Mnr; Sr];
else;
    Y=(Hd-repmat(Mnd,mH,1)) ./ repmat(Sd,mH,1);
    Y=(Y.*repmat(Sr,mH,1))+1;
    S=[nanmean(Y);  nanstd(Y)];
end;



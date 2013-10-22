function stats=validat1(y,X,b,I2,yrscal,yrsval,kopt)
% validat1:  validation of a stepwise forward regression model
% s=validat1(y,X,b,I2,yrscal,kopt);
% Last revised 5-16-01
%
% You estimated a regression model by forward stepwise entry, say using stepr1.m.  How well does the model validate 
% at each step of the stepwise procedure?  You might use this information to select a simpler model, or possibly just
% to quantify the likely predictive value of the the regression on independent data
%
%*** INPUT
%
% y (my x 2)r predictand, with year as column 1
% X (mX x nX)r   time series  matrix of potential predictors, year as col 1
% b (1 x ?)r   regression coefficients, constant first, then coefficients on predictors, in same order as in I2
% I2 (1 x ?)i  col index to X indicating which variables (cols) of X entered at each step of regression
%   Note that must strip year col off X before applying the column index I2
% yrscal (1 x 2)i  start and end year of calibration period
% yrsval (1 x 2)i  start and end year of validation period, or [] if method is cross-validation
% kopt (1 x 1)i  options
%    kopt(1) method of validation
%       ==1 cross-validation
%       ==2 split sample
%
%
%*** OUTPUT
%
% stats-- structure variable with statistics (? steps in the stepwise regression)
%
%   .R2 (? x 1)r   calibration R-squared at each step of the entry;  last value should match your regression R-squared
%   .MSE (? x 3)r   mean square error of calibration (col 1), validation (col 2), and validation for null recon (col 3)
%   .RMSE (? x 2)r  root mean square error of calibration and validation 
%   .RE (? x 1)r  reduction of error statistic
% 


% Pull calib data
Lc = y(:,1)>=yrscal(1) & y(:,1)<=yrscal(2);
yc = y(Lc,2);
yryc = y(Lc,1);

yrX = X(:,1);
X(:,1)=[];
Lc =yrX>=yrscal(1) & yrX<=yrscal(2);
Xc = X(Lc,:);
yrXc = yrX(Lc);

if ~all (yrXc == yryc);
    error('Mismatch in calibration years');
end;

% Compute calib-period mean of predictand
ycmean = mean(yc);
ncal = length(yryc);

nfull = length(I2); % number of predictors in final or fill model
nsteps = length(I2);


alpha=0.95;  % unneeded

if kopt==1; % if crossvalidate
    Lx = crospull(ncal,0,0); % logical matrix of years for cr
    
    for m =1:nsteps; % loop over steps -- crossvalidate at each
        
        I3=I2(1:m); % pointer to predictors
        
        % Cross check on calibration using all data
        XXc = [ones(ncal,1) Xc(:,I3)];
        [B,BINT,Rcal,RINT,STATS] = REGRESS(yc,XXc,alpha);
        R2cal = STATS(1);
        
        
        
        yhatcv = repmat(NaN,ncal,1); % to hold crossvalid estimates
        ecv = repmat(NaN,ncal,1); % to hold crossvalid errors
        for n = 1:ncal; % loop over years
            Lthis = Lx(:,n);
            Xthis = [ones(ncal-1,1)  Xc(Lthis,I3)];
            ythis = yc(Lthis);
            Lout = ~Lthis;
            yout = yc(Lout);
            xout= [1 Xc(Lout,I3)]; % predictor data for deleted residual
            [b,BINT,R,RINT,s] = REGRESS(ythis,Xthis,alpha);
            yhout = xout*b;
            yhatcv(n)=yhout;
            ecv(n)=yout-yhout;
        end;
        %   .R2 (? x 1)r   calibration R-squared at each step of the entry;  last value should match your regression R-squared
        %   .MSE (? x 3)r   mean square error of calibration (col 1), validation (col 2), and validation for null recon (col 3)
        %   .RMSE (? x 2)r  root mean square error of calibration and validation 
        %   .RE (? x 1)r  reduction of error statistic
        % 
        stats.R2{m}=R2cal;
        MSEc = (sum(Rcal .^2))/(ncal-m-1) ; % residual mean square
        MSEv = (sum(ecv .^2))/ncal; % mean square error of cross-v
        MSEn = sum((yc-ycmean) .^2)/ncal;  % mean square error for null recon
        stats.MSE{m}=[MSEc MSEv MSEn];
        stats.RMSE{m}=[sqrt(MSEc) sqrt(MSEv) sqrt(MSEn)];
        stats.RE{m}=  1 - (MSEv/MSEn);
        
        
    end; % for m=1:nseps
    
else;  %if kopt==1; % if crossvalidate

end;% if kopt==1; % if crossvalidate

    
    

        

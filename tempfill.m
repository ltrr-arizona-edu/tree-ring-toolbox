function est=tempfill(Y,X,overlap,yr,m)

% Estimation of missing monthly temperature value by extrapolated
% departure-from-mean method (Conrad and Pollack 1950)


%****   INPUT ARGUMENTS  ******
%
% Y   (m1 x 13) monthly temp array for station needing estimation
%		First col is year.
% X   (m2 x 13) monthly temp array for station to be used as estimator
% overlap (1 x 2)  first and last year of overlapping period to
%		be used for computing "long-term" means.
% yr	year of value to be estimated
% m   month of value to be estimated
%

%*******  OUTPUT ARGUMENT *******
%
% est  (1 x 1) the estimated value for year yr and month mo
%



I=find(X(:,1)==yr);  % index to year for estimation


% Form 0-1 vectors pointing to overlap period in Y,X 
LY = Y(:,1) >= overlap(1)  & Y(:,1) <= overlap(2);
LX = X(:,1) >= overlap(1)  & X(:,1) <= overlap(2);

% Pull out subset of overlap years for the specific month of est
% Compute long-term mean monthly temperatures for each
%		station for overlap period.
% Assign same departure-from-mean at predictor station to missing
%		value at predictand station.

Y1 = Y(LY,m+1);
X1 = X(LX,m+1);

%**** New code to allow using periods with missing data to be used in 
%**** computing long-term means 

LY1=Y1>=9998;
LX1=X1>=9998;;
LL=LY1 | LX1;
Y1(LL)=[];
X1(LL)=[];

%******  end of new code


xm = mean(X1);  % long-term mean, or mean for overlap period
ym = mean(Y1);
dep=X(I,m+1) - xm;  % Departure from mean at station x in target year
est = ym+dep;

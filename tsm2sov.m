function [v,YRS]=tsm2sov(X,kopt,yrs)
%
%  [v,YRS]=tsm2sov(X,kopt)
%
% time series matrix to strungout vector
%
% Meko 2-19-97
%
%********************* IN *******************************88
%
% X (mX x nX)r time series matrix, mX years, nX cols. Maybe a years 
%		column, depending on kopt(1)
% kopt (1 x 1)i   options
%		kopt(1) is there a leading year column in X
%			==1 no
%			==2 yes
% yrs (1 x 2)i   start, end year of X. Required only if kopt(1)==1
%
%********************** OUT ********************************
%
% v (mv x 1)r the strung out vector 
% YRS (nsers x 3)i  start year, end year, and row index into v of
%		the time series data for each series stored in v
%
%************************* OTHER FUNCTIONS CALLED 
%
% intnan.m  -- internal nan in a vector?
%
%********************** NOTES *************************************
%
% So far, requires that no "internal" nans in any series, just leading or
% trailing nans


% Size
[mX,nX]=size(X);




% Get year vector for X, and if needed strip years column off X
if kopt(1)==1; % no year column
	if exist('yrs')~=1;
		error('No yrs variable passed as input arg')
	else
		[ mtemp,ntemp]=size(yrs);
		if mtemp~=1 | ntemp~=2,
			error('yrs must be 1 x 2');
		end
		yr1 = (yrs(1):yrs(2))';
	end
else; % X has a years column
	yr1=X(:,1);
	d1=diff(yr1);
	if ~all(d1==1);
		error('Supposed years column in X not increment by 1');
	end
	X(:,1)=[];
end

[mX,nsers]=size(X);  % nsers is number of time series

% Check that no internal NaNs in any series
for n=1:nsers;
	x=X(:,n);
	mcheck=intnan(x);
	if mcheck==1;
		error(['Internal NaN in series ' int2str(n)]);
	end
end


% Size storage matrix for years and pointer
a=NaN;
YRS=a(ones(nsers,1),ones(3,1));
mlength = mX * nsers;
v=a(ones(mlength,1),:);

% Pointer to  NaNs in X
L1=isnan(X);


% Loop over series, building the sov and the years reference matrix
isp=0;
for n=1:nsers;
	x=X(:,n);
	L2=L1(:,n);
	x(L2)=[]; % strip of leading and trailing NaNs
	yrx=yr1;
	yrx(L2)=[]; % year vector for stripped series
	nyears = length(yrx);

	igo=isp+1;
	isp=igo+nyears-1;
	v(igo:isp)=x;  
	YRS(n,:)=[yrx(1) yrx(nyears) igo];
end;

% Strip off extra NaNs from v
L3=isnan(v);
if any(L3);
	v(L3)=[];
end






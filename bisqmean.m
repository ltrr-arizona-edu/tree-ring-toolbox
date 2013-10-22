function [ymn,varyh,df,w,ybar,se]=bisqmean(y)
%
% Biweight mean for a vector of numbers.
% D Meko 2-18-95
%
% Source:  Mosteller and Tukey (1977, p. 205, p 351-352)
%		Cook and  Kairiukstis (1990, p. 125-126)
%
%
%****************  INPUT *************************
%
% y (? x 1)r  vector of data -- say, indices for ? cores in a year
% 
%
%********************  OUTPUT ************************
%
% ymn (1 x 1)r  biweight mean
% s (1 x 1)r   asymptotic standard dev of biweight mean - p. 208, 
%		third eqn from top of page
% w (? x 1)r  final weights on values in y
% ybar (1 x 1)r arithmetic mean corresponding to ymn
% se (1 x 1)r  standard error of ybar
%
%********************  NOTES *********************
%
% ybar and v1 just included for debugging purposes to double check
% on closeness of ybar to ymn, v1 to v
%
%*******************************************************************


sens = 0.001;  % hard coded theshold of sensitivity for stopping iterat
nits = 10;  % max number of allowed iterations

[n,ny]=size(y);
if ny > 1;
	error('y should be a vector')
end
if n<6, error('Should use median for n<6'), end


ww = 1/n; % weight for even average
ybar = mean(y); % arith mean
se= sqrt(var(y)/n); % standard error of mean

nz=0;
ymn = ybar; % initial biweight mean as arith mean

for i = 1, nits;  % iterate max of nits times
	ymnold = ymn;  % store old value of mean
	e = y-ymn; % deviations from mean
	S = median(abs(e));  % median abs deviation
	u = e / (6*S);  % scaled deviations

	w = (1 - u.^2).^2;  % compute weights
	L1 = u>=1; % flag huge errors
	L1s = sum(L1);
	if L1s>0
		nz=0;
		nz= nz(ones(L1s,1),:);
		w(L1)=nz;  % set weights on those obs to zero
	end
	w = w / sum(w); % adjust weights to sum to 1.0
	
	ymn = sum(w .* y); % compute biweight mean


	% Variance of estimate of biweight mean
	ui= e / (9*S);
	L2 = ui>1;
	ui(L2)=[];
	z =y(~L2);
	nz = length(z);
	nom1 = (z - ymn) .^2;
	nom2 = (1-ui .^2) .^4;
	nom = sum(nom1 .* nom2);

	den1 = sum((1-ui .^2) .* (1-5*ui .^2));
	den2 = -1 + sum ((1-ui .^2) .* (1-5*ui .^2));
	varyh = nom / (den1*den2);  % variance of biweight mean 
	% last eqn, p. 208
	df = 0.7 * (nz -1); % degrees of freedom 


   % if little change in mean, exit loop
	if abs (ymn - ymnold) < sens
		break
	end
end


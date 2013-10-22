function k=subyears(pdall,pdsht)

% Compute subperiods of long term equal in length to short period
% Use: Have reconstruction with a calibration period.  Want to
%    divide the reconstruction into subperiods equal in length to
%    the calibration period -- say, to compare frequency of drought
%    in modern vs previous periods

nsht=pdsht(2) - pdsht(1) + 1;
nall=pdall(2) - pdall(1) + 1;


nbefore=pdsht(1)-pdall(1);  % number of years before modern period
nc=floor(nbefore/nsht);  % number of full periods before modern period
left = rem (nbefore,nsht); %  number of years in odd first subperiod

yrbeg= (pdall+left:nsht:pdsht(1))';
yrend=  yrbeg+nsht-1;

k = [yrbeg yrend];


if left~=0
	k1=  [pdall(1)  pdall(1)+left-1];
	k = [k1 ; k];
end

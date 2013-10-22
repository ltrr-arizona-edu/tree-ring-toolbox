function x=trailnan(x);
% trailnan: remove any trailing NaN's from a vector
% x=trailnan(x);
% Last revised: Meko 2-28-97
%
% trailnan is a utility function used by rwlinp.m and other functions. Most users
% will not need to call trailnan in their own code 
%
%*** INPUT
%
% x (1 x ?)r or (? x 1)r vector, usually a time series
%
%*** OUTPUT
%
% x (? x 1)r  the same vector, but with any trailing NaN's removed
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- none
%*** NOTES -- none

[mx,nx]=size(x);
if mx==1 & nx>1; % row vector -- to  cv
	x=x';
	mx=nx;
	nx=1;
elseif mx>1 & nx==1; % col vector -- ok
else
	error('x must be cv or rv');
end


L1=isnan(x);
sum1=sum(L1);
if sum1==0,
	return
else; % at least 1 NaN
	if ~isnan(x(mx));  % last value in x is not NaN, so no trailing NaN
		return
	else; % have at least one trailing nan
		fgood = find(~isnan(x)); % row index to valid data
		fhigh=max(fgood); % index of last valid value
		iout=(fhigh+1):mx; % want to remove these trailing elements
		x(iout)=[];
	end
end

function [igo,isp]=consec(x)
% [igo,isp]=consec(x) 
%
% Longest consecutive sequence of valid (not NaN) data in a 
% vector
%
% D Meko   5-5-96
%
%********************** IN ARGS *******************************
%
% x (mx x nx)r vector -- either mx or nx should be 1
%
%
%*********************** OUT ARGS ***************************
%
% igo -- subscript of first  value in the good sequence
% isp -- subscript of the last value in the sequence
%
%************************* NOTE *****************************
%
% If all values in x are NaN, returns  igo=0 and isp=0


[mx,nx]=size(x);
if mx~=1 & nx ~=1, 
	error('x must be a row or column vector')
end

% Convert to rv if needed
if nx==1;
	x=x';
end

yr = (1:(length(x)+2)); % dummy "years" vector

if all(isnan(x))
	igo=0;
	isp=0
elseif all(~isnan(x))
	igo=1;
	isp=length(x);
else
	L1=~isnan(x);
	L2 = [0 L1 0];
	d=diff(L2);
	i1 = find(d==-1);
	i2 = find(d==1);
	i3 = i1-i2;  % lengths of NaN periods
	[m,im]=max(i3);
	isp = yr(i1(im))-1;
	igo = isp - m +1;
	n = isp - igo+1;
end



	


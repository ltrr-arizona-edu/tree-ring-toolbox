function [R,phi]=polar(x,y)
%
% [R,phi]=polar(x,y)
%
% Convert real and imaginargy parts of a series of complex numbers to
% their magnitudes and phases.  
%
%  D Meko 7-20-95
%
%
%********************  IN ARGS *******************************
%
% x (mx x 1)r real part of series
% y (yx x 1)r imaginary part of series
%
%
%********************* OUT ARGS *********************************
%
% R (mx x 1)r magnitude
% phi(mx x 1) phase -- possible range is -pi to pi
%
%
%******************** METHOD *********************************
%
% Serves same purpose as subroutine "polar" in Bloomfield (1976, p. 150)
%
%

[mx,nx]=size(x);
[my,ny]=size(y);
if nx~=1 & ny~=1,
	error('x and y must be col vectors')
end
if mx ~=my,
	error('x and y must be same length')
end

R = sqrt(x .^2 + y .^2);
phi = atan2(y,x);

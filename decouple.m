function z= decouple(y,x)
% dcouple:  decouple two variables; remove linear dependence on a second variable
% CALL z = decouple(x,y);
%
% Meko 10-8-98
%
%********* IN 
%
% y(? x 1)r  time series of first variable
% x(? x 1)r  time series of second variable
%
%
%********** OUT
%
% z(? X 1)r  residual from linear regression of y on x  
%
%******** NOTE
%
% Example.  Have temperature sereis y and precip series.  In correlating precip and
% temperature with tree rings, want to remove confusing influence of correlation
% between pcp and tmp.  Set y to tmp, x to pcp.  Get z, that part of the temperature
% unpredictable from linear relationship with pcp.  In Panofsky notation, this is
% T.P

%--------  CHECK INPUT

[mx,nx]=size(x);
[my,ny]=size(y);

if nx~=1 | ny~=1;
   error('x and y must be column vectors');
end
if mx ~= my;
   error('x and y must be same length');
end


if any(isnan(x)) | any(isnan(y));
   errror('x and y not allowed to have NaNs');
end

%---------  want constant term in model
xmean=mean(x);
ymean=mean(y);

x=[ones(mx,1) x];

[b,bint,z,rint,stats]=regress(y,x,.01);

% restore mean to residuals
z = z + ymean;

ymean;
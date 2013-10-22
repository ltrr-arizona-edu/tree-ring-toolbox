function qplot(x,f);
% qplot: quantile plot
% qplot;
% Last revised 8-24-00
%
% Quantile plot of data
%
%*** INPUT
%
% x (mx x 1)r  data
% f (mf x 1)r <optional> f-values for quantiles
%
%*** OUTPUT
%
% No output args.
% 
% Plot is made of f-th uantiles as function of f.  Vertical lines at median and
% first and third quartiles.
%
%*** REFERENCES
%
% Cleaveland, W. S., 1993.  Visualizing Data.  Hobart Press, Summit, New Jersey, 360 p.
%   [P 16-20] 
%
%*** UW FILES CALLED
%
% qtile.m
%
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES 
%
% See notes to qtile.m for info on allowable f



if any(isnan(x));
   error('x must contain no NaNs');
end;

[mx,nx]=size(x);
if nx~=1; 
   error('x must be col vector');
end;
n=mx;


if nargin==1;
   f=[];
   [q,fout]=qtile(x);
   f=fout;
else;
   [q,fout]=qtile(x,f);
   if length(fout)~=length(f);
      error('Length of fout not eq to length of x');
   end;
   
end;

[mf,nf]=size(f);
str1=sprintf('n_x=%5.0f',n);
str2=sprintf('n_q=%5.0f',mf);


hp1=plot(f,q,f,q,'.');
xlabel('f');
ylabel('q(f)');
title(['Quantile Plot; (' str1 ';  ' str2 ')']);
grid;
set(gca,'XLim',[0 1],...
   'XTick',[0 0.25 0.50 0.75]);

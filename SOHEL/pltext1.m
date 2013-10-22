function pltext(x,y,txt)
%
% USAGE : pltext(x,y,txt)
%   Adds text txt to a plot at decimal fraction x from left side and
%   decimal fraction y from bottom.
%
%
% INPUTS
%-------
% x (1 x 1)	x coordinate as a fraction of x-axis range (0-1).
% y (1 x 1)	y coordinate as a fraction of y-axis range (0-1).
% txt		String variable to be placed on the figure window.
%
%
% NO OUTPUTS
%
% 
% NO USER WRITTEN FUNCTIONS NEEDED
%_________________________________________________________________

v=axis;  % gets [xmin xmax ymin ymax] of current axes

x1=v(1) + x * (v(2)-v(1));   % x plotting point
y1=v(3) + y * (v(4)-v(3));   % y plotting point

text(x1,y1,txt,'Fontsize',8)
%axis;

% end of function
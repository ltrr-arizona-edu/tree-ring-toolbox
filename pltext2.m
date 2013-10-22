function pltext2(x,y,yup,fontsize,txt)
% plot text at specified position; paired function with setylim2.m
% D Meko 10-17-94
%
%
%******************************** IN ARGS *******
%
% x (1 x 1)r   start text this far from the far left of the x axis,
%	where this far is given as decimal fraction of total 
%	range covered by x axis
% y (1 x 1)r   likewise for y-axis value at which text is to begin;
%	but decimal fraction here refers to only that part of
%	y-axis range from yup to the upper  y limit.
% yup (1 x 1)r   lower threshold of y-axis of reserved text area
%	(originally get with a call to setylim1.m from calling fctn)
% fontsize (1 x 1)i size of font
% txt (1 x ?)s   text
%
%**********************  OUT ARGS -- NONE   ************
%
%*********************  USER-WRITTEN FUNCTIONS NEEDED -- NONE ****
%
%*********  NOTES ******************************************
%
% Example
%
% pltext2(.1,.5,.8,8,'This')
%
% Says plot the text string 'This' 0.1 of the way from the left 
% xlimit and half (0.5) of the way between a y-value of 0.8 and the
% upper y limit.  Use font size of 8.
%
% The text is centered vertically on the calculated y position


a = axis; % get current axis limits
xpos = a(1) + x*(a(2)-a(1));
ypos = yup + y*(a(4)-yup);
text(xpos,ypos,txt,'FontSize',fontsize) 
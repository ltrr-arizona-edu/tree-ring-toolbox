%
% USAGE : CLSGRF
%   This is a script .m file containing a set of 
%   commands for the uicontrol callback function 
%   to close all figure windows
%___________________________________________________________

close(1);

while gcf ~=1,
  close(gcf);
end

close(1);

% End of file

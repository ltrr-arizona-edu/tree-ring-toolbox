function [cL,cB,cW,cH]=figsize(fwidth,fheight);
% figsize: get acceptable computer-specific figure position coordinates
% [cL,cB,cW,cH]=figsize(fwidth,fheigth);
% Last Revised 2-9-04
%
% Get acceptable computer-specific figure position coordinates.  User
% specifies desired figure width and height as decimal fraction of the
% screen width and height.
%
%****   INPUT 
%
% fwidth (1 x 1)r   figure width as decimal fraction of screen width
% fheight (1 x 1)r  figure height as decimal fraction of screen height
%
%*** OUTPUT
%
% cL (1 x 1)i  left position as pixel
% cW (1 x 1)i  width in pixels
% cB (1 x 1)i  bottom position as pixel
% cH (1 x 1)i  height in pixelsesult -- structure of results
%
%*** REFERENCES --- none
%
%*** TOOLBOXES NEEDED -- none
%
%*** UW FUNCTIONS CALLED -- none
%
%*** NOTES
% 
% In calling program, to resize figure window use:
% set(gcf,'Position',[cL cB cW cH])

Ssz=get(0,'screensize'); % Get screen size


if fwidth<0.1 | fwidth>0.9;
    error('Width as fraction of screen width must be between 0.1 and 0.9');
end;
if fheight<0.1 | fheight>0.9;
    error('Height as fraction of screen height must be between 0.1 and 0.9');
end;


cW=floor(fwidth*(Ssz(3)-Ssz(1)));
cH =floor(fheight*(Ssz(4)-Ssz(2)));
cL = floor(0.5* (Ssz(3)-cW));
cB = floor(0.5* (Ssz(4)-cH));


if (cL+cW)>=Ssz(3);
    error('Screen width cannot accomodate specified position');
end;
if (cB+cH)>=Ssz(4);
    error('Screen height cannot accomodate specified position');
end;


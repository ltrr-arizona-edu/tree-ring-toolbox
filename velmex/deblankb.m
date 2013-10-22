function s=deblankb(s);
% deblankb:  deblank a string on both left and right
% s=deblankb(s);
% Last revised 10-11-01
%
% Deblank a string on both left and right
%
%*** INPUT 
%
% s (1 x ?)s  string
%
%
%*** OUTPUT
%
% s (1 x ?)s deblanked string
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%*** NOTES


if ~isstr(s);
    error('s must be string');
end;


s=fliplr(deblank(fliplr(deblank(s))));
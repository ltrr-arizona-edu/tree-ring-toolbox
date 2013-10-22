function [y,xinc]=shiftup1(x);
% [y,xinc]=shiftup1(x);
% shiftup1: shift a non-negative time serie by adding smallest non-zero value 
% Last revised 1-1-01
%
% Utility function for shifting a time series prior to operations that would be invalid on zero values.
%
%*** INPUT
%
% x (1 x mX)r time series, must be non-negative, but may have zero values
%
%*** OUTPUT
%
% y (mY x 1)r  optionally shifted time series
% xinc (1 x 1)r amount of shift applied to x to produce y (see notes)
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Input series is checked for zero values.  If any zeros, the series is shifted by incrementing by the
% smallest nonzero value


[mX,nX]=size(x);
if nX~=1;
    error('x must be col vector');
end;

% Check x
L=  (x<0);
if any(L);
    i1=find(L);
    E=[i1 x(i1)];
    disp(' i   x(i)');
    disp(E);
    error('x must be non-negative; displayed values invalid');
end;

L=x==0;
if all(L);
    error('all values of x are zero');
end;

if any(L);
    x1=x;
    x1(L)=[]; % remove zero values
    xinc=min(x1);
else;
    xinc=0;
end;

y = x + xinc;
    
    
    
    
    
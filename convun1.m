function Y=convun1(X,dtype,unitin,unitout,multin,multout)
% convun1: convert units of climatic data
% CALL: Y=convun1(dtype,unitin,unitout,multin,multout);
%
% Meko 5-9-97
%
%******************** IN 
%
% X (mX x nX)r input climatic data
% dtype (1 x ?)s  data type. Example: PPT
% unitin (1 x ?)s units of input.  Example: in (inches)
% unitout (1 x ?)s units of output. Example: in (inches)
% multin (1 x ?)r  multiplier of input.  Example: 0.01 (input is 100ths of inches)
% multout ( x ?)r  multiplier of output.  Examaple: 1.00 (output is inches)
%
%********************* OUT 
%
% Y (mY x nY)r  output climate data, in converted units
%
%******************** NOTES 
%
% convun1.m converts by the equation y = (a*(x+b)+c)/d
% settings for a, b, c, d are determined from the input parameters.  For example,
% if the input x is in hundredths of inches and you want output in inches, 
% would have  b=0,c=0,a=0.01, d=1.0
%
% Bookkeeping on years is assumed to be handled in the calling program.  In other
% words, convun1.m simply converts a matrix of numbers from one type of units to
% another.  In a typical application, X might have 12 cols, representing the
% months jan to dec.  But X could also be a row vector, column vector, or even
% a scalar.  X is allowed to have NaNs, which will be NaNs in Y.
%
% Conversions implemented so far:
%
% dtype: ppt, tmp
% unitin,unitout: mm,in   --- F,C


%----------------- Check Input

if nargin~=6;
   error('convun1 requires 6 input args');
end

if ~isnumeric(X) | ~ischar(dtype) | ~ischar(unitin) | ~ischar(unitout);
   error('X must be numeric, and dtype, unitin, and unitout must be char');
end
if ~isnumeric(multin) | ~isnumeric(multout);
   error('multin and multout must be numeric')
end

X=X*multin;


switch dtype
case {'PPT','ppt'};
   switch unitin;
   case {'in','IN'};
      switch unitout;
      case {'in','IN'};
         c=0; b=0; d=1.0;
         a=1.0
         
      case {'mm','MM'};
         c=0; b=0; d=1.0;
         a=25.4;
      otherwise
         error([unitout ' is invalid unitout']);
      end; % of switch unitout
   case {'mm','MM'};
      switch unitout;
      case {'mm','MM'};
         c=0; b=0; d=1.0;
         a=1.0;
      case {'in','MM'};
         c=0; b=0; d=1.0;
         a=(1/25.4);
      otherwise
         error([unitout ' is invalid unitout']);
      end
   otherwise
      error([unitin ' is invalid unitin']);
   end
   
case {'tmp','TMP'};
     switch unitin;
   case {'F','oF'};
      switch unitout;
      case {'F','oF'};
         c=0; b=0; d=1.0;
         a=1.0;
      case {'C','oC'};
         c=0; b=-32; d=9.0;
         a=5.0;
      otherwise
         error([unitout ' is invalid unitout']);
      end; % of switch unitout
   case {'C','oC'};
      switch unitout;
      case {'C','oC'};
         c=0; b=0; d=1.0;
         a=1.0;
      case {'F','oF'};
         c=32; b=0; d=1.0;
         a=9/5;
      otherwise
         error([unitout ' is invalid unitout']);
      end
   otherwise
      error([unitin ' is invalid unitin']);
   end
   
     
otherwise
   error([dtype ' is invalid dtype']);
end


% Convert
Y = (a*(X+b)+c)/d

Y=Y/multout;
   
      
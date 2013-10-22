function [xsegm,YM,yrY,fl] = xycolw(xseg,Y,yY)
%
% USAGE : [xsegm,YM,yrY,fl] = xycolw(xseg,Y,yY)
%   Given two time series x and y, builds matrices for use in 
%   statistical functions relating a segment of x to all same 
%   length segments of y.
%
%
% INPUTS
%-------
% xseg(mx x 1)	Time series xseg.
% Y (ny x 1)	Time series Y.
% yY (ny x 1)	Year vector for Y time series.
%
%
% OUTPUTS
%--------
% xsegm (seglen x ?)	Duped columns of segment of x defined by
%                       segend, seglen.
% YM (seglen x ?)	Offset (by one year) segments of y equal
%                       in length to seglen.
% yrY (? x 1)		Ending years of each series (col. of) in Y
% fl (1 x 1)		Error flag. If there is error, fl=-1.
%
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% UDISP.M 	An interactive display function.
%________________________________________________________________________
  
  seglen=length(xseg);
  if (length(Y) < seglen);
    udisp('Length of Y should be greater than or equal to length of xseg.');
    fl=-1;
    return;
  end

  cols = length(Y) - seglen + 1;
  xsegm = xseg(:,ones(cols,1));

  A = (1:seglen)';
  A = A(:,ones(1,cols));
  B = (0:cols-1);
  B = B(ones(seglen,1),:);
  Y1 = A + B;
  Y1= Y(Y1);
  YM = zeros(seglen,cols);
  YM(:) = Y1;
  yrY = (yY(1)+seglen-1 : yY(length(Y)))';                 

% End of function 

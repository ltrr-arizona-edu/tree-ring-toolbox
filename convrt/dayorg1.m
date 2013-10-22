function [y,ycode,yrnew,mnew]=dayorg1(c,kopt,yrold,mold)
% dayorg1: subfunction; organize a month of daily Eischeid P data
% [y,ycode,yrnew,mnew]=dayorg1(c,kopt,yrold,mold);
% Last revised 12-5-00
%
% Called by day2fle1 to organize a line of input file. 
%
%*** IN
%
% c ( 1 x 
% kopt
%   (1)==1 first record of new site
%      ==2 another record of ongoing site
% yrold (1 x 1)i  previous year (e.g. 1978);  [] if kopt(1)==1
% mold (1 x 1)i  previous month (e.g., 1==Jan); [] if kopt(1)==1
%
%
%*** OUT
%
% y (1 x 29,30 or 31)r  daily P for each day of month
% ycode ( 1 x 29,30 or 31)i code for daily values
% yrnew (1 x 1)i  current year
% mnew (1 x 1)i current month
%
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES USED -- none


% Get current year,month
yrnew=str2num(c(1:4));
mnew=str2num(c(5:6));
nday=str2num(c(8:9));

% Year and month consistent with previous?
if kopt(1)==2; % not the first record for this stn
   if mnew==1; % if current month Jan
      if mold~=12;
         error(['mnew is 1, but mold is not 12']);
      end;
      if yrnew~=yrold+1;
         error(['yrold, yrnew = ' str2num([yrold yrnew])]);
      end;
   elseif mnew<=12 & mnew>0 ;
      if mnew~=mold+1 | yrnew~=yrold;
         error(['yrold,  mold, yrnew, mnew = ' str2num([yrold mold yrnew mnew])]);
      else;
         % all OK;
      end;
   end;
end;



% Compute start and end index of data
if nday<28 | nday>31;
   error('nday must be between 28 and 31');
end;
ion=13; % start char index
ioff=12+nday*7; % end index

% Pull data
c1=c(ion:ioff);
c2=(reshape(c1,7,nday))';


% Store output
ycode=str2num(c2(:,7));  % estimation code in position 7
y = str2num(c2(:,1:6));

if size(y,1)==28; 
    if (mnew==2 & rem(yrnew,4)==0) & rem(yrnew,100)==0  & rem(yrnew,400)~=0;
        y(29)=NaN;
        ycode(29)=NaN;
    elseif ~(mnew==2 & rem(yrnew,4)~=0);
        error('28 days in month & Feb of a leap yr');
    else;
        y(29)=NaN;
        
        ycode(29)=NaN;
    end;
end;

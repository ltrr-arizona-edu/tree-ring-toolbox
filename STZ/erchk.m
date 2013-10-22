function [erflg,xind1,xind2]=erchk(yrv,first,last)
  
% Get the row indices of specified beginning and ending years of a segment
% of a time series whose full length covers the years in vector yrv
% And make sure the specified endpoints are consistent with yrv
%
% Called by cfit1.m


erflg=0; % error flag

% Check first year of segment
if first<=yrv(1); % Default to beginning year of series
  xind1=1;
elseif first>=yrv(length(yrv)); % Start year pointed to is equal or
	% later than last year in time series
  erflg=-1;
  return;
else
  xind1=find(yrv==round(first)); % accept nearest year as beginning year
end

% Check last year of segment
if last>=yrv(length(yrv)), 
  xind2=length(yrv);
elseif last<=yrv(1),
  erflg=-1;
  return;
else
  xind2=find(yrv==round(last));
end 

% End of file

function [mnx,stdx,n] = basic1(yrs,IX)
 
% DESCRIPTION : [mn,std,n] = basic1(yrs,IX);
% Calculates and returns the mean and standard deviation of 
% the core indices
% 
% INPUTS  :  yrs (ns x 3) - Year matrix
%	     IX (? x 1)   - Strung out vector containing indices
%
% OUTPUTS :  mnx (ns x 1)  - Mean of each core IDs
%	     stdx (ns x 1) - Standard deviation for each core IDs
%	     n (ns x 1)    - Number of years for each core ID without
%			     considering NaNs.
%_____________________________________________________________________

% Initialize 
[ns,dms]=size(yrs);
dum1=NaN;
dum2=dum1(ones(ns,1),:); % vector of NaN

mnx=dum2;
stdx=dum2;
n=dum2;

% Cull out the indices for each core sequentially
for scid=1:ns,
 %yrv=(yrs(scid,1):yrs(scid,2))';
 xv=IX(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));

 lzn=isnan(xv);
 zn=xv;
 zn(lzn)=[];
 if ~isempty(zn),
  mnx(scid)=mean(zn);
  stdx(scid)=std(zn);
  n(scid)=length(zn);
 end
end


% End of file

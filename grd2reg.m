function [P,RA,R,RM,RS,JP3,wt2]=grd2reg(IR2,CG,Q,PM,PS,JP2)
% [P,RA,R,RM,RS,JP3,wt2]=grd2reg(IR2,CG,Q,PM,PS,JP2)
% 
% gridpoint monthly climate data to regional; also gridpoint
% stdzd anomalies to gridpoint monthly time series in climatic
% units; and in a non-regional mode, just converts the gridpoint
% stdzd anomalies to original climate units at the gridpoints.
% 
% Meko 2-6-97; Last revised 4-27-00
%
%********************** IN ARGS ********************************
%
% IR2(mreg x nreg)i index to rows of CG telling which gridpoints
%		are in which region.  Each col is for a region. Cols
%		padded with zero as needed to matrix. In non-regional mode,
%		IR2==[];
% CG (mG x 2)r  long/lat mapping coordinates, mG gridpoints
% Q (mQ x nQ)r multi-station matrix of monthly standardized
%		anomaly time series for gridpoints.  mQ years by
%		nQ columns.  First 12 cols of Q are Jan-Dec data for
%		gridpoint 1, etc.  
% JP2 (mG x 12)i  index to cols of Q telling which cols are
%		Jan data for all gridpoints, Feb data, etc. For example,
%		to get all Jan data, take the columns pointed to by the
%		first col of JP2.
% PM (mG x nG)r long-term means, Jan-Dec, for gridpoints
% PS (mG x nG)r long-term std devs
%
%******************** OUT ARGS *******************************
%
% P (m4 x n4)r like Q, but converted from stdzd anomalies to
%		original climate units
% RA (m5 x n5)r regional stdzd anomalies computed by weighting
%		gridpoint anomalies proportional to cos of latitude. Col
%		size of RA is nreg x 12, where nreg is number of regions
% R (m5 x n5)r like RA, but in original climate units. R is 
%	computed by (1) getting the latitude weighted long-term 
%	means and standard deviations, and (2) using those to convert
%  RA.
% RM (nreg x 12)r long term means for regions. 
% RS (nreg x 12)r long term means for regions. 
% JP3 (nreg x 2)i index to cols of R,RA telling first, last
%	 cols to grab data for region 1, 2, ... nreg
% wt2(mG x 1)r  weights on gridpoints in CG
%**************** NOTES **************************************
%
% grd2reg.m is a dedicated function specifically to be called by
% regcli2.m.  Refer to regcli2.m for descriptions of the input.


% Get number of gridpoints
[mG,nG]=size(CG);
npts=mG;

% Size
[mQ,nQ]=size(Q);
[m2,n2]=size(JP2);
if m2~=mG | n2~=12,
	error('Pointer mtx JP2 should be same row-size as gridpoint xy matrix')
end

% Allocate
P=repmat(NaN,mQ,nQ);

if npts==1; % One gridpoint or target site; no regional series to be computed
   RA=[];   R=[]; RM=[]; RS=[]; wt2=[]; nR=[];
   if ~isempty(IR2);
      error('Non-regional mode requires that input arg IR2 be []');
   end;
else; % One-gridpoint-at-a-time mode; no regionalization
   [mreg,nreg]=size(IR2); % nreg is number of regions
   nR=nreg*12; % number of cols in R, RA
   RA=repmat(NaN,mQ,nR);
   R=repmat(NaN,mQ,nR);
   RM=repmat(NaN,nreg,12);
   RS=repmat(NaN,nreg,12);
   wt2=repmat(NaN,mG,1); % cv of weights on each gridpoint 
end;

% Build pointer JP3 to allow retrieving individual regional means and stdevs from  RM, RS
if npts==1; 
   JP3=[];
else;
   j1=(0:(nreg-1))';
   jgo=j1*12+1;
	jsp=jgo+11;
   JP3=[jgo jsp]; % 2-cols; col 1 is start col index for region, col 2 is end col index
end;

% If more than 1 gridpoint, loop over regions
if npts > 1;
   for n = 1:nreg;
      
      % get vector of gridpoint pointer
      ir2=IR2(:,n);  
      % remove any trailing zeros
      ir2(ir2==0)=[];
      num2=length(ir2); % number of gridpoints in this region
      
      % get vector of latitudes for the gridpoints; get weights
      glat=CG(ir2,2); 
      % convert to radians, and take cos
      wt=cos(glat * pi/180);
      % make wts sum to 1.0
      wt = (wt /sum(wt))'; % cv of weights
      wt2(ir2)=wt;
      
      % Pull the set of long-term means for region, std devs
      XM=PM(ir2,:);
      XS=PS(ir2,:);
      
      % Weight the gridpoint long-term means and std devs 
      RM(n,:)=wt*XM;
      RS(n,:)=wt*XS;
      
      % Weight the time series of monthly stdzd gridpoint anomalies
      % Loop over months
      for k=1:12;
         % Compute col in R to store this region/month's data
         rcol=k + (n-1)*12;
         
         % Grab data for this month of year, all gridpoints
         E = Q(:,JP2(:,k));
         % Cull data for this region's points
         E = E(:,ir2);
         
         % Weight the gridpoint stdzd anomalies to regional
         e=E*wt'; % regional values this region/month, all years
         RA(:,rcol)=e;
         
         % Calc the corresponding data converted to original units
         rm=RM(n,k); %  long-term mean this region/month
         rs=RS(n,k); % long-term std dev
         f = e*rs + rm;  % transform to climatic units
         R(:,rcol)=f;
      end;
   end;
end;



% Convert the gridpoint stdzd anomaly series to climatic units

% String out PM as row vector
PM1=  PM';
pm1=(PM1(:))';
% dupe rv to matrix, same size as Q
PM1=pm1(ones(mQ,1),:);

% String out PS as row vector
PS1=  PS';
ps1=(PS1(:))';
% dupe rv to matrix, same size as Q
PS1=ps1(ones(mQ,1),:);

% Multiply elements by std dev, add mean
P = (Q .* PS1) +PM1;


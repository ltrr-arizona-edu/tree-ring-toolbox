function [P,n,s]=runs(x,yr,xc,ksign);
% runs: runs analysis -- years of runs, run length, run sum
% CALL: [P,n,s]=runs(x,yr,xc,ksign);
%
% Meko 4-23-97
%
%********************* IN **************************************
%
% x (mx x 1)r time series.  No NaN allowed
% yr (mx x 1)i year vector for x
%  
% xc (1 x 1)r threshold critical value of time series
% ksign (1 x 1)i   ==1 : event is x>xc
%                  ==-1: event is x<xc
%
%
%***************************  OUT **********************
%
% P (nruns x 2)i start and end year for each run
% n (nruns x 1)i length of each run
% s (nruns x 1)i run sums

%******************Check inputs

a=NaN;
[mx,nx]=size(x);
if nx~=1;
   error('x must be cv')
end
[mtemp,ntemp]=size(yr);
if mtemp ~=mx | ntemp~=1;
   error('yr must be cv, same size as x');
end

[mtemp,ntemp]=size(xc);
if mtemp~=1 | ntemp~=1;
   error('xc must be scalar');
end
[mtemp,ntemp]=size(ksign);
if (mtemp~=1 | ntemp~=1) | (ksign~=-1 & ksign~=1);
   error('ksign must be -1 or 1');
end


% Get start, end years of runs
Ls=runstart(x,xc,ksign);
Le=runend(x,xc,ksign);
Ls=logical(Ls);
Le=logical(Le);
nruns=length(Le);
if length(Ls) ~=length(Le)
   error('Ls and Le must be same length');
end
P=a(ones(nruns,1),ones(2,1)); % initialize

% Get run length
n=runleng(x,Ls,Le);

% Compute periods for runs
P=[yr(Ls) yr(Le)];

% Get run sum
s=runsum(x,Ls,Le,xc,ksign);

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


function L=runstart(x,xc,ksign);
% runstart: logical pointer to start years of runs in a time series
% CALL: L=runstart(x,xc,ksign);
%
% Meko 4-23-97
%
%********************* IN **************************************
%
% x (mx x 1)r time series.  No NaN allowed
%  
% xc (1 x 1)r threshold critical value of time series
% ksign (1 x 1)i   ==1 : event is x>xc
%                  ==-1: event is x<xc
%
%
%***************************  OUT **********************
%
% L (mx x 1 )L  1 if value starts a run, 0 if not

%******************Check inputs

[mx,nx]=size(x);
if nx~=1;
   error('x must be cv')
end
[mtemp,ntemp]=size(xc);
if mtemp~=1 | ntemp~=1;
   error('xc must be scalar');
end
[mtemp,ntemp]=size(ksign);
if (mtemp~=1 | ntemp~=1) | (ksign~=-1 & ksign~=1);
   error('ksign must be -1 or 1');
end


% Get first mx-1 and last mx-1 values of x
x1=x(1:(mx-1)); % first segment
x2=x(2:mx); % last segment; first value is at time 2, last value at time mx

if ksign==-1;
   L1=x1>=xc & x2<xc; % x2 value starts a run
   % Handle start year
   if x(1)<xc;
      L=[1 ; L1];
   else
      L=[0 ; L1];
   end
else
   L1=x1<=xc & x2>xc;
   % Handle start year
   if x(1)>xc;
      L=[1 ; L1];
   else
      L=[0 ; L1];
   end
end


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

function L=runend(x,xc,ksign);
% runend: logical pointer to end years of runs in a time series
% CALL: L=runend(x,xc,ksign);
%
% Meko 4-23-97
%
%********************* IN **************************************
%
% x (mx x 1)r time series.  No NaN allowed
%  
% xc (1 x 1)r threshold critical value of time series
% ksign (1 x 1)i   ==1 : event is x>xc
%                  ==-1: event is x<xc
%
%
%***************************  OUT **********************
%
% L (mx x 1 )L  1 if value ends a run, 0 if not

%******************Check inputs

[mx,nx]=size(x);
if nx~=1;
   error('x must be cv')
end
[mtemp,ntemp]=size(xc);
if mtemp~=1 | ntemp~=1;
   error('xc must be scalar');
end
[mtemp,ntemp]=size(ksign);
if (mtemp~=1 | ntemp~=1) | (ksign~=-1 & ksign~=1);
   error('ksign must be -1 or 1');
end


% Get first mx-1 and last mx-1 values of x
x1=x(1:(mx-1)); % first segment
x2=x(2:mx); % last segment; first value is at time 2, last value at time mx

if ksign==-1;
   L1=x1<xc & x2>=xc;
   % Handle end year
   if x(mx)<xc;
      L=[L1; 1];
   else
      L=[L1; 0];
   end
else
   L1=x1>xc & x2<=xc;
   % Handle end year
   if x(mx)>xc;
      L=[L1;  1];
   else
      L=[L1;  0];
   end
end


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

function m=runleng(x,Ls,Le);
% runleng: run lengths for all runs in a time series
% CALL: L=runleng(x,Ls,Le);
%
% Meko 4-23-97
%
%********************* IN **************************************
%
% x (mx x 1)r time series.  No NaN allowed
% Ls (mx x 1)L  1 if corresponding year in x is start of a run
%               0 otherwise
% Le (mx x 1)L  1 if corresponding year in x is end of a run
%               0 otherwise
%
%
%***************************  OUT **********************
%
% m (nruns x 1 )i  length of each run in the time series

%******************Check inputs

[mx,nx]=size(x);
if nx~=1;
   error('x must be cv')
end

% Get start and end index for runs
igo = find(Ls);
isp = find(Le);
nruns = length(igo);
if length(igo) ~= length(isp);
   error('row size of igo and isp should be same');
end

m  = isp - igo +1;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

function s=runsum(x,Ls,Le,xc,ksign);
% runsum: run sum for all runs in a time series
% CALL: L=runsum(x,Ls,Le);
%
% Meko 4-23-97
%
%********************* IN **************************************
%
% x (mx x 1)r time series.  No NaN allowed
% Ls (mx x 1)L  1 if corresponding year in x is start of a run
%               0 otherwise
% Le (mx x 1)L  1 if corresponding year in x is end of a run
%               0 otherwise
% xc (1 x 1)r threshold critical value of time series
% ksign (1 x 1)i   ==1 : event is x>xc
%                  ==-1: event is x<xc
%
%
%
%***************************  OUT **********************
%
% s (nruns x 1 )i  run sums

%******************Check inputs
a=NaN;


[mx,nx]=size(x);
if nx~=1;
   error('x must be cv')
end

% Get start and end index for runs
igo = find(Ls);
isp = find(Le);
nruns = length(igo); % number of runs
if length(igo) ~= length(isp);
   error('row size of igo and isp should be same');
end

s=a(ones(nruns,1),:);



% Make logical pointer to observations in each run
L1=zeros(mx,nruns);
for n=1:nruns;
   j=igo(n):isp(n);
   L1(j,n)=1;
end
L1=logical(L1);

% Compute run sums
if ksign==1; % events are values greater than xc
   d=x-xc;
else
   d=xc-x;
end
for n=1:nruns;
   LL1=L1(:,n);
   s(n)=sum(d(LL1));
end

function runsnoi1
% runsnoi1: runs analysis via runs.m on reconstruction with added noise
% CALL: runsnoi1;
%
% Written for Jordan pcp reconstruction analysis.  Needed to get runs info on
% the actual data as well as the reconstruction with noise as intermediate
% step before runsbox1.m 
%
% Meko 12-12-97
%
% ************** IN ***********************
%
% No input arguments.
%
% User points to two files:
%
% 1 - <rn*.mat> with the following series stored from a reconstruction
%
%   yrec (? x 2)r reconstructed predictand , year in col 1
%   yact (? x 2)r actual predictand; typically for the calibration period
%      but could extend beyond that
%   YRS (2 x 2)i start and end year of 
%      row 1: calibration period used for reconstruction
%      row 2: desired summary period (typically, reconstruction + calibration period)
%   ycrit(1 x 1)r  critical value defining an event in yrec and yact
%   ksign(1 x 1)i  tells whether an event is value greater or less than ycrit
%      ==1 greater
%      ==-1 less
%   fnsim (1 x ?)s  path/filename to file holding the noise series, typically
%      generated previously by recnois1.m
% 2 - <runnout?.m>  output file to hold the following
%   d1 (? x 1)i run lengths for calibration period actual predictand data
%   d2 (? x 1)i run lengths for calibration period reconstruction
%   d4 (? x 1)i run lengths for noiseless reconstruction -- including the
%		part in the calibration period as well as the earlier part
%   D(? x nsims)i  run lengths for each of the simulations -- which are recons
%		with added noise; these source time series come from file fnsim
%	 s1,s2,s4,S:  run sums corresponding to the above run sums
%
%********************  OUT ****************
%
% No output arguments
%
% Output is the .mat file described above 
%


%--------------  GET INPUT DATA--------------------------

[file1,path1]=uigetfile('rn*.mat','Input data for runsnoi1.m');
pf1=[path1 file1];
eval(['load ' pf1]);
icheck=[exist('yrec') exist('yact') exist('YRS') exist('ycrit')...
      exist('ksign') exist('fnsim')];
if  ~all(icheck==1);
	error(['All required variables not in ' pf1]);
end


%----------------- GET FILENAME FOR OUTPUT STORAGE FILE

[file4,path4]=uiputfile('*.mat','Output storage file');
pf4=[path4 file4];



%-----------------  ALLOW USER TO OVERRIDE NUMBER OF GROUPS FOR OUTPUT FILE
titdlg = 'User can override default n of groups for runbox1.m input';
prompt = {'Enter the number of groups'};
def = {20};
lineno = 1;
ngroup=inputdlg(prompt,titdlg,lineno,def);
ngroup=ngroup{1};

% -------------- CHECK INPUT DATA

% yact
[mtemp,ntemp]=size(yact);
if ntemp~=2;
	error('yact should be 2-col tsm');
end
nyr1 = mtemp;
yract=yact(:,1);

% yrec
[mtemp,ntemp]=size(yrec);
if ntemp~=2;
	error('yrec should be 2-col tsm');
end
nyrrec = mtemp;
yrrec = yrec(:,1);

% YRS
[mtemp,ntemp]=size(YRS);
if ntemp~=2 | mtemp~=2;
	error('YRS must be 2 x 2');
end

% YRS must be consistent with years in yact and yrec
if YRS(1,1) < min(yract);
	error('YRS(1,1) too early for yact');
end
if YRS(1,2) > max(yract) | YRS(1,2)>max(yrrec);
	error('YRS(1,2) too recent for either yact or yrec');
end
if YRS(2,1)<min(yrrec);
	error('YRS(2,1) before first value in yrrec');
end
if YRS(2,2)> max(yrrec);
	error('YRS(2,2) too late for yrrec');
end


% ycrit
[mtemp,ntemp]=size(ycrit);
if mtemp~=1 | ntemp~=1;
	error('ycrit must be scalar');
end

% ksign
[mtemp,ntemp]=size(ksign);
if mtemp~=1 | ntemp~=1;
	error('ksign must be scalar');
end
if ksign~=-1;
   error('Not yet ready to handle wet events');
end



% fnsim;
if ~ischar(fnsim);
	error('fnsim must be string');
end



%------------- FORM POINTERS


L1c = yract>= YRS(1,1) & yract<=YRS(1,2); % to calib period in yact --variable 1
L2c = yrrec>= YRS(1,1) & yrrec<=YRS(1,2); % to calib period in yrec -- variable 2
L2r = yrrec>=YRS(2,1) & yrrec<= YRS(2,2); % to recon+noise period in yrec -- vbl 4


% ------------  PULL TIME SERIES AND YEAR VECTORS

y1 = yact(L1c,2);  % actual data for calib period
yr1 = yact(L1c,1);
nyr1 = length(yr1);
y2 = yrec(L2c,2); % recon data for calib period
yr2 = yrec(L2c,1);
nyr2 = length(yr2);
y4 = yrec(L2r,2); % recon+noise analysis period
yr4 = yrec(L2r,1);
nyr4 = length(yr4);

% means and variances of departures
defic = -1.0 * (y1-ycrit);
Ltemp = defic>0;
defic(~Ltemp)=NaN;
mn1 = nanmean(defic);
var1 = (nanstd(defic)) .^2;

defic = -1.0 * (y2-ycrit);
Ltemp = defic>0;
defic(~Ltemp)=NaN;
mn2 = nanmean(defic);
var2 = (nanstd(defic)) .^2;

defic = -1.0 * (y4-ycrit);
Ltemp = defic>0;
defic(~Ltemp)=NaN;
mn4 = nanmean(defic);
var4 = (nanstd(defic)) .^2;


%--------------  GET THE RECON+NOISE SERIES

eval(['load ' fnsim]);
if ~exist('Y')==1;
	error(['Y does not exist in ' fnsim]);
end

[nyr4,nsim]=size(Y);
if nyr4~=length(yr4);
	error(['row size of Y in ' fnsim ' not equal to length of yr4']);
end

% Compute mean and variance of annual deficits
if ksign==-1;
   defic = -1.0 * (Y-repmat(ycrit,nyr4,nsim));
   Ltemp = defic>0;
   defic(~Ltemp)=NaN;
   mnns = nanmean(defic);
   varns =  (nanstd(defic)) .^2;
else
   error('ksign==1 not yet valid');
   
end
   

%-------------- GENERATE THE RUNS INFO -----------------------

% Actual data, calibration period
[P,d1,s1]=runs(y1,yr1,ycrit,ksign);
n1 = length(d1); % number of runs


% Recon data, calibration period
[P,d2,s2]=runs(y2,yr2,ycrit,ksign);
n2 = length(d2); % number of runs


% Recon data, full period
[P,d4,s4]=runs(y4,yr4,ycrit,ksign);
n4 = length(d4); % number of runs


%--------------- ALLOCATE FOR  STORING NOISY RUN RESULTS

D = repmat(NaN,nyr4,nsim); % for run lengths
S = repmat(NaN,nyr4,nsim); % for run sums

%----------------- RUNS ANALYSIS OF NOISY RECONS

for m = 1:nsim; % loop over simulations
	z = Y(:,m); % the time series
	[P,d,s]=runs(z,yr4,ycrit,ksign);
	num1 = length(d); % number of runs
	if num1>0;
		D(1:num1,m)=d;
		S(1:num1,m)=s;
	else
		error('Hey, no runs in any of simulated series');
	end
end
	

%------------------  STORE RESULTS

set1 = [' D d1 d2 d4 S s1 s2 s4 nyr1 nyr2 nyr4 ngroup'];
set2 = [' mn1 mn2 mn4 mnns ']; % means of annual deficits
set3 = [' var1 var2 var4 varns ']; % variance of annual deficits
eval(['save ' pf4 set1 set2 set3]);



 



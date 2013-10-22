function surrgt1
% surrgt1: replace missing calib-pd data in predictor matrix with surrogate series
% CALL: surrgt1
%
% Meko 4-21-97
%
%**********************  IN ************************************************
%
% 1. User prompted .mat file name.  File has time series matrix X and
%    year vector yr.  Note that this file was probably built with sov2tsm2.m,
%    which also built an ascii file with 1 line per tree giving information
%    the cols of X.  Usually will want to have this reference ascii file handy.
% 2. User prompted .mat file with additional required input
%
%    yrscal (1 x 2)i first and last year of calibration period
%    mlap (1 x 1)i  minumum number of years overlap for accepting a correlation 
%      coefficient between series.  These coefficients used for selecting
%      surrogate series.
%    J1 (? x 2)i <empty> specified surrogate series; empty if none.
%      Col 1 points to a col of X and denotes a series
%      Col 2 points to another col of X and denotes a series to use as surrogate
%    I1 (? x ?)i grouping of the series in X into subgroups
%      Each row corresponds to a separate group. For example, matrix X might
%      have trees from the rincons, catalinas and galiuros (3 groups).
%   
%      Each col gives pointer to cols of X indicating the series to include n
%      in the group.  For example:
%
%       [1 2 3 4 5 NaN NaN;
%        6 7 8 9 NaN NaN NaN; etc
%
%      says that cols 1,2,3,4, and 5 in X belong to group 1, and
%      cols 6,7, 8 and 9 to group 2
%
%***************************** OUT *********************************************
%
% 1. User prompted .mat filename.  File like input file #1, but with the revised
%    time series matrix as X.  Year vector is still yr.
% 2. User prompted .txt filename.  Information on the surrogates. One line per
%    series needing a surrogate. Example
%
%    1(0.65/48)-3   
%    2(0.55/59)-6-3
%      etc
%
%    In this example, line 1 says that that series 3 was used as the surrogate for 
%    series 1, and that the correlation between 1 and 3 based on 48 yr of data 
%    is 0.65.    The (1) indicates that series 1 is an "order 1" series.
%    Line 2 says that series 2 , of order 2,  is daisy chained to calibration series 3
%    by way of series 6, and that the correlation between 2 and 6 is 0.55, based
%    on 59 yr of data. See subfunction bldstr1.m for more on the ascii file.
%
%    "Order" is defined as follows:
%
%     0 series overlaps entire calibration period
%     1 daisy chained directly to an order-0 series
%     2 daisy chained to an order-0 series by way of an order-1 series
%     3 daisy chained to an order-0 series by way of an order-2 series, which connects
%        to an order-1 series, which connects to an order-0 series
%     etc
%     99 series has an 'assigned' surrogate, meaning that correlations not used
%        to identify surrogate
%
%     "overlap" means a common period of data of at least mlap years

a=NaN;
Jbig=[];
Dbig=[];
ordbig=[];

% Get .mat file with X and yr
[file1,path1]=uigetfile('*.mat','.mat file with input tsm X and year vector yr');
pf1=[path1 file1];
eval(['load ' pf1]);
if ~exist('X')==1 | ~exist('yr')==1;
   error('X or yr not in .mat input file')
end


% Check size of X, and that row size consistent with  yr
[mX,nX]=size(X);
[mtemp,ntemp]=size(yr);
if mX ~= mtemp | ntemp~=1;
   error('X must have same rowsize as yr, and yr must be a cv');
end

% Rename number of years and variables in X
nsers = nX;
nyrs = mX;

% Allocate for the revised matrix Y, and for cv of order of series
Y=X; % initialize output matrix as the input matrix
ordbig=a(ones(nsers,1),:);

%******************* Get file with the additional input

[file2,path2]=uigetfile('*.mat','.mat file with additional input');
pf2=[path2 file2];
eval(['load ' pf2]);
if ~exist('J1')==1 | ~exist('mlap')==1 | ~exist('yrscal')==1 | ~exist('I1')==1;
   error('J1, I1, mlap or yrscal not in .mat input file')
end

% Check J1
if isempty(J1);
   mJ1=0;
   nJ1=0;
else
   [mJ1,nJ1]=size(J1);
   if nJ1~=2;
      error('Col size of J1 must be 2')
   end
end



% Compute number of groups; build YRS
ngroups = size(I1,1);
YRS=[min(yr) max(yr); yrscal];

% Initialize cv to hold number of problem series in each group
NNP=a(ones(ngroups,1),:);


%********************* MAPPING OF SERIES TO SUBGROUP NUMBERING 
A=a(:,ones(nsers,1));
% Find max number of series in any group
L1=~isnan(I1);
nmax =  max(sum(L1'));

%Loop over first to maximum number nmax
for n = 1:nmax;
   % get a col of I1 and strip off NaNs
   i1=I1(:,n);
   i1=i1(~isnan(i1)); % this is a rv
   % dupe the number n into a rv same size
   nvec = n(:,ones(length(i1),1));
   % Fill A slots
   A(i1)=nvec;
end
% check that all elements mapped
if any(isnan(A));
   error('A still has NaN elements')
end
% The rv A now maps the columns of I to columns of group matrices

 
%***************************  LOOP OVER GROUPS ***********************
np=0; % running cumulative counter of series needing surrogates
for n =1:ngroups;
	nnp=0; % running counter of series needing suurogates in each group
	disp(['Group ' int2str(n)]);
   i1 = I1(n,:); % which series are  in the group; pointer to col of X
   i1=i1(~isnan(i1));
   num1=length(i1);  % number of series in group
   
   % Build prespecified surrogate matrix JS
   if isempty(J1);
      JS=[];
   else
      II1 = i1(ones(mJ1,1),:);
      j1=J1(:,1);
      JJ1=j1(:,ones(num1,1));
      L1=II1 == JJ1;
      L2 = (any(L1'))';
      JS=J1(L2,:);
      JS=[A(JS(:,1)) A(JS(:,2))];
   end
   
   
   % Pull subset of the data corresponding to this group
   X1=X(:,i1);
   
   % Fill in the calib period segment of submatrix with surrogate data
   [Y1,J,D,ord,nump]=surrgt2 (X1,YRS,JS,mlap,rmin,1);
	nnp=nnp+nump; % increment counter
   np=np+nump; % increment cumulative counter of number of problem series
   
   % Fill series order into ordbig
   ordbig(i1)=ord;
   
   % Put the revised data in correct cols of Y
   Y(:,i1)=Y1;
   
   % Map col pointers to X1 in J to cols of X & augment Jbig
   if ~isempty(J);
      [mJ,nJ]=size(J);
      Jnew = a(ones(mJ,1),ones(nJ,1));
      for k = 1:2; % loop over cols 1 and 2 of J
         jj=J(:,k); % a cv
         L1=~isnan(jj);
         if sum(L1)>0;
            jjj=jj(L1);
            jjj=(i1(jjj))';
            Jnew(L1,k)=jjj;
         end
      end
      Jbig=[Jbig; Jnew];
   end
   
   % Map pointer in col 1 of D to X instead of X1
   % Note that col 1 of D should be same as col 1 of J in call to surrgt2.m
   if ~isempty(D);
      D(:,1) = Jnew(:,1);
      Dbig=[Dbig; D];
   end
NNP(n)=nnp;
end; % of loop over groups
      

%******************************* OUTPUT THE .MAT DATA STORAGE  ************
[file3,path3]=uiputfile('*.mat','.mat outfile to hold filled in tsm');
pf3=[path3 file3];
X=Y; % overwrite X with Y.  No further need for original X, and want name of 
  % stored output matrix to be X
eval(['save ' pf3 ' X yr ']);

%***************************** OUTPUT THE ASCII FILE ************************

% Check that number of problem series matches row size of Jbig and Dbig
if np~=size(Jbig,1) | np~=size(Dbig,1),
   error('np not equal to row size of Jbig and Dbig');
end

% Open file for output
[file4,path4]=uiputfile('*.txt','.txt outfile of daisy chain info');
pf4=[path4 file4];
fid4=fopen(pf4,'w');

% Loop over problem series
nnp=0;
j=1;

for n =1:np;
	ntarget=NNP(j);
   strout1=bldstr1(ordbig',Jbig,Dbig,n);
   fprintf(fid4,'%s',strout1);
	nnp=nnp+1;
	if nnp==ntarget;
		fprintf(fid4,'%s\n','****************************************');
		nnp=0;
		j=j+1;
	end
	
end; % for n=1:np
fclose (fid4);

   
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


function [Y,J,D,ord,nump]=surrgt2 (X,YRS,JS,mlap,rmin,kopts)
% surrgt2.m: surrogate tree-ring indices  
% CALL: [Y,J,D,ord,nump]=surrgt2 (X,YRS,JS,mlap,rmin,kopts)
%
%************************* IN ******************************************
%
% X(mX x nX)r time series matrix, mX years by nX variables.  Any col
%	may contain NaN's, but the non-NaN part must be contiguous: i.e., 
%	no imbedded internal NaNs in any series
% YRS(2 x 2)i start and end year of
%		row1: matrix X
%		row2: calibration period -- period that want to fill with
%			data for all series
% JS (mJS x 2)i column index to X specifying preordained surrogate series
%		for some series. 
%		col 1: pointer to col of X, specifying first series
%		col 2: pointer to col of X, specifying surrogate series
%		 Example of need.  Have core indices for various
%		trees at a site.  Know that tree7 was right next to tree8 and looked
%		exactly the same.  Tree7 covers 1300-1850, and tree8 covers 1600-1997.
%		Decided that rather than rely on correlations, will merely specify
%		on logical reasoning that tree8 is a good surrogate to fill in the
%		1910-1997 calibration period data for tree7.
%		JS==[] means no preordained surrogates
% mlap (1 x 1)i minimum acceptable overlap for computing correlation coef
% rmin (1 x 1)r  minimum acceptable correlation coef
% kopts(1 x 1)i options
%		kopts(1) Fast mode or debugging mode
%			==1 fast
%			==2 debugging
%
%*********************** OUT ********************************************
%
% Y (mY x nY)r time series matrix like X, but with surrogate data filled in
%		for previous NaNs in calibration period
% J (mJ x 2)i summary information on diasy chained surrogate series
%		col 1: col in X holding the key series
%		col 2: col in X holding the chosen immediate surrogate series
%	
%		Note that J is returned as [] if all series are order-0
%
%		Note terminology for order of the daisy chaining:
%			0=is a calib series (overlaps entire calibration period)
%			1=series whose surrogate is order-0
%			2=series whose direct surrogate is order-1
%			3=series ... order-2
%			-1= prespecified surrogate
%
% D (mD x 3)r  correlation coef and sample size (# yr) between
%		key series and its immediate series in the daisy chain.  Row size of
%		D is same as row size of J; each row corresponds to a ps.  Columns:
%			1 - key series pointer
%			2 -- correlation coef
%			3 -- sample size for the correlation
% ord (? x 1)i order of each of the series in X. See above for order def
% nump (1 x 1)i number of problem series (non calibration series)
%
%  Note: ultimate goal in the calling program surrgt1.m is to
%  be able to give ascii record of the daisy chaining. For each problem 
%  series (row), will want to list the seq number, order, 
%  and seq nos and orders of all series back to the order-0 series 
%  used as surrogate. Example:
%
%		3-6-2  r=.65 (120) 1921-1990
%
%	indicates that series #3, of order 2, was daisy chained back to 0-order
%	series #2 by way of 1-order series #6;  and that the correlation between
%	series #3 and #6 based on 120 years of data is 0.65. D is retured as []
%  if no surrogates were needed 
%
%***************** APPROACH **************************************
%
% cs=calibration series: overlaps calib period completely
% ps= problem series: others
%
% surrgt2.m establishes heirarchy of daisy chain
% starting with the correlation matrix.  
% Assume time series matrix X, corresp years info on X period and
% calibration period; and optionally, specified surrogates
% 
% Make a copy of the original tsm X -- say the copy is Y
%
% Convert X to a sov, using tsm2sov.m
%
% Compute correlations of each series with all other series with r4sov.m, giving
% main results in R
% 
% Calc col indices of calibration period series (cs series), defined as those
% with data in X  covering entire calibration period.  Store this and col indices
% of other series, defined as problem series (ps series)
%
%
% Build col pointers to cs and ps: jc1, jp1
% Initialize and allocate for string matrix D
% Initilize variable pointers jc2 to jc1, jp2 to jp1
% 
% Update Y,jc2,jp2,ord,order1, order2, order3, for any pre-specified surrogates
% If at least mlap overlap years in the key and surrogate, use call to scale1.m.
% If fewer than mlap overlap years, or if no overlap, use call to scale2.m, whic
% scales and shifts relative to global means and standard deviations for each
% series.
%
% If all ps series handled, return.
%
% While loop until all ps handled
%	Check for highest correlation between any ps and cs or revised ps
%		The cs series is the surrogate
%		Set the order of the revised ps to 1 if the 'done' series order 0,
%			or to 2 if the 'done' series order 1, etc.
%		Store the pointer to the surrogate, the correl coef, and sample size
%		Update Y by replacing key series segment with its surrogate
%		Update pointers jc2,jp2, 
%		Update ord
%		Update J by appending row with [nkey jsur]
%		Update D by appending row with [nkey, r, sample size]
%		Exit while if all ps handled
%		Error mssg if some ps left and no correlation avail (insuff overlap
%			with any other series daisy chained to cs
%	End of while	
%	
%
%
% Note:  surrgt1.m is the calling function.  In surrgt1.m have tsm with
% series from possibly many sites. Some of these series will be grouped 
% together.  For example, study  is of san pedro basin, with core indices
% from galiuros, rincons, santa ritas.  Separate groupings and corresponding
% calls are made to surrgt2.m for each grouping.  Surrgt1.m handles bookeeping
% and revising of the master tsm. This bookeeping includes cross-references of
% column indices, and printing out of a text file documenting the daisy
% chain
%   


%**********************  USER FUNCTIONS CALLED ***********************
%
% scale1.m -- shift and scale times series to desired mean and variance
% surrgt3a.m -- returns revised series given original series and surrogate
% series
%	(calls scale1.m)
% surrgt3b.m -- like surrgt3a.m, but for pairs without overlap of mlap obs
% tsm2sov.m -- time series matrix to strung out vector

% Initialize variables
a=NaN;
D=[];
J=[];

%--------------------------- Check inputs

if nargin~=6;
	error(['Needs 6 input args, has ' int2str(nargin)]);
end

if kopts(1)~=1 & kopts(1)~=2
	error('kopts(1) must be 1 or 2');
end

[nyrs,nsers]=size(X);
vseq=(1:nsers); % rv of col numbers of X

[mYRS,nYRS]=size(YRS);
if mYRS~=2 | nYRS~=2
	error('YRS must be 2 x 2')
end
yrc = (YRS(2,1):YRS(2,2))'; % year vector for calibration period
yr=(YRS(1,1):YRS(1,2))'; % year vector for X
if max(yrc)>max(yr) | min(yrc)<min(yr);
	error('Specified calib period not in period covered by X');
end



if isempty(JS);
	% no specified surrogates
else
	if size(JS,2) ~=2;
		error('Col size of JS must be 2');
	end
	if any(any(JS>nsers));
		error('An element of JS exceed max possible seq number of series');
	end
	if size(JS,1)>nsers;
		error('Row size of JS greater than nsers');
	end
end


%--------------  Make some year vectors

yr=(YRS(1,1):YRS(1,2))'; % full period of X
yrc=(YRS(2,1):YRS(2,2))'; % calibration period



%*************** MAKE A COPY OF DATA MATRIX X

Y=X;



%**************** MAKE A SOV VERSION OF X
% r4sov.m uses sov as input 
[v,vyrs]=tsm2sov(X,1,YRS(1,:));


%  store the calibration period start and end years
p=YRS(2,:);

%**********************  COMPUTE THE CORRELATIONS
%
% Call r4sov.m to get the correlations of series with all other series
% In the call to r4sov.m, v is the sov, vyrs is the corresponding matrix of
% start year, end year, and start row index of series in v;  mlap is the 
% minimum acceptable overlap to accept a correlation coefficient; p is a 
% 1 x2 vector specifying the "calibration" period, and "1" is an option setting
% specifying that the correlations are to be computed on the full overlap period
% rather than just the overlap period within the calibration period

R = r4sov(v,vyrs,mlap,p,1);


%--------------- Build col pointers to cs and ps: jc1, jp1
%
% cs means "calibration series" -- those overlapping entire calibration period
% ps means "problem series" -- the other series
% jc1 is pointer to cols of X specifying which series are cs series
% jp1 is pointer to cols of X specifying which series are ps series

% Initialize row vectors of pointers
jc1=[];
jp1=[];
order1=[];
order2=[];
order3=[];
ord=a(:,ones(nsers,1));  % initialize order number for all series in X
num1=0; num2=0; num3=0;

% Get matrix of calib-period segment of X
L1=yr>=YRS(2,1) & yr<=YRS(2,2);
X1=X(L1,:);

% Make column pointer to order-0 series in X
L2=all(~isnan(X1));  % rv
jc1=vseq(L2); % pointer to calib series
% compute number of order-0 series
if ~isempty(jc1); % have at least one cs
   numc=sum(L2); % number of cs
   ord(jc1)=0;  % set elements of ord for these series to "zero" order
else; % have no cs
	numc=0;
	% All surrogates should have been specified
	if size(JS,1)~=nsers;
		error('No order-0 series and all surrogates not specified');
	end
end


% Make column pointer to all other series (problem series ) in X
jp1=vseq(~L2); % pointer to problem series
% compute number of ps series
if ~isempty(jp1); % have at least one ps
   nump=sum(~L2);
else
	nump=0; % no problem series
	return
end



% Make logical pointer to any specified series
jj3=logical(zeros(1,nsers));
if isempty(JS);
   % no action needed, j3 is the pointer
else
   numjs = size(JS,1);
   jj3(JS(:,1))=1;
end


% Allocate storage matrices for problem series info

J = a(ones(nump,1),ones(2,1)); % see opening comments
D = a(ones(nump,1),ones(3,1)); % see opening comments



% Initialize running identifier of remaining problem and order-0 series
jp2=jp1;
jc2=jc1;


% Initialize logical pointers
jj1=logical(zeros(1,nsers)); % for order-0 series
jj2=jj1;  % for problem series
jj1(jc2)=1;
jj2(jp2)=1;
jj3=jj3; % built previously

jj2a=jj2 & ~jj3;  % points to problem series excluding specified series


%******************  WHILE LOOP TILL SURROGATES BUILT
%  FOR ALL EXCEPT SPECIFIED SERIES

k1=1; % while control
kount1=0; % counter for problem series - for filling D and J
while k1
   if ~any(jj2a); % problem series all taken care of
      k1==0;
      break
   else
      [ikey,isurr,rsurr,nsurr,errflg]=highr2(R,jj1,jj2a,jj3,mlap,rmin);
      if errflg==0;
         kount1=kount1+1;
         D(kount1,:)=[ikey rsurr nsurr];
         J(kount1,:)=[ikey isurr];
         ord(ikey)=ord(isurr)+1;
         jj2a(ikey)=0;
         jj2(ikey)=0;
         jj1(ikey)=1; 
         
         % Update time series matrix
         y=X(:,ikey);
         x=Y(:,isurr);
         % Get the revised form of y
         yrev = surrgt3a(x,y,yr,yrc);
         % Substitute the revised y into a column of Y
         Y(:,ikey)=yrev;
         
      end; % if errflg==0
   end; % if ~any(jj2a)
end; % while k1


%************************************************************
% Problem series now all done, except for any with specified
% surrogates. Now handle those, if any.
if ~any(jj3); % no specified surrogates
   if any(jj2); 
      error('Still have unfinished series');
   end
else; % some specified series to handle
   ord(jj3)=99;
   % Loop over j3 series
   numj3=sum(jj3);
   for n=1:numj3;
      i1=JS(n,1);
      i2=JS(n,2);
      y=X(:,i1);
      x=Y(:,i2);

      kount1=kount1+1;
      J(kount1,:)=[i1 i2];
      
      L1=R(:,3)==i1 & R(:,4)==i2;
      L2=R(:,3)==i2 & R(:,4)==i1;
      
      if any(L1 | L2);  % Sufficient overlap for this pair
         
         % Get the correlation
         if sum(L1)==1 & sum(L2)==0;
            D(kount1,:)=[i1 R(L1,1) R(L1,2)];
         elseif sum(L2)==1 & sum(L1)==0;
            D(kount1,:)=[i1 R(L2,1) R(L2,2)];
         else
            error('More than one row of R has this pair');
         end
         
         
         % Update time series matrix
         yrev = surrgt3a(x,y,yr,yrc);
         % Substitute the revised y into a column of Y
         Y(:,i1)=yrev;
         
      else; % No correlations for this pair
         
         D(kount1,:)=[i1 NaN NaN];
         
         % Update time series matrix
         yrev = surrgt3b(x,y,yr,yrc);
         % Substitute the revised y into a column of Y
         Y(:,i1)=yrev;
      end;   % any(L1 | L2) 
   end;  %  n=1:numj3;
end; % ~any(j3); % no specified surrogates

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

function yrev=surrgt3a(x,y,yr,yrc);
%
% Compute surrogate data for part of one time series, given overlap
% data for that series with another.  Called by surrgt2.m
%
% Meko 2-21-97
%
% 
%******************IN
%
% x(mx x 1)r one time series
% y(my x 1)r another series, needing estimate
% yr (m1 x 1)i year vector for x, y
% yrc (m2 x 1)i year vector for calibration period.  Any data in this
%   period that is NaN in y will be replaced with estimated data
%		based on stdzd anomalies using stats from overlap of x and y
%
%******************* OUT
%
% yrev (my x 1)r revised y
%
%**************** CALLS
%
% scale1.m
%
%******************** NOTES ******************************
%
% Uses method described in scale1.m
% Dedicated function for surrgt2.m
%


[mx,nx]=size(x);
[my,ny]=size(y);
if mx ~=my;
	error('x and y must be same length');
end
if nx~=1 | ny~=1;
	error('x and y must be col vects');
end


% Get overlap segments of x,y
L1=~isnan(x) & ~isnan(y);
x1=x(L1);
y1=y(L1);


% Make pointer to part of calib period needing data for y
yrgo=yrc(1);
yrsp=yrc(length(yrc));
L2=yr>=yrgo & yr<=yrsp & isnan(y);

% Get that data for x, and check that no NaN in this data
x2=x(L2);
if any(isnan(x2)); 
	error('NaN in the predictor series x2 in calib period');
end

% Call scale1.m to get estimated data for that period in y
y2=scale1(x1,y1,x2);


%  Substitute the estimates in y
y(L2)=y2;
yrev = y;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


function yrev = surrgt3b(x,y,yr,yrc);
% surrgt3b:  surrogate series for assigned surrogates; no overlap version
% CALL: yrev = surrgt3b(x,y,yr,yrc);
%
% Meko 4-20-97
%
%************************ IN 
%
% x (mx x 1)r  predictor series (surrogate series)
% y (my x 1)r  predictand series (series needing data)
% yr (myr x 1)i   year vector, applies to x and y
% yrc (myrc x 1)i  year vector for calibtation period
%
%
%*************************** OUT
%
% yrev (mx x 1)r revised for of y, with NaN data in calibration
%  period replaced with 'standardized anomaly' estimates
%
%*************************** NOTES **********************
%
% x has complete data for calibration period.  y has some
% or all data NaN in calibration period.  Objective is to
% replace the calibration period NaN of y with estimates from
% x.  Estimates are such that standardized anomalies from 
% the long-term mean are the same for the missing-data period
% for series x and y.


[mx,nx]=size(x);
[my,ny]=size(y);

if nx~=1 | ny~=1;
   error('col size of x, y should be 1');
end
if mx~=mx;
   error('x and y should have same row size');
end
if length(yr) ~= length(x)
   error('yr must be same length as x, y');
end

% Initialize yrev as y
yrev=y;

% Year settings
yrgo=yr(1);
yrsp=yr(mx);

yrcgo=yrc(1);
yrcsp = yrc(length(yrc));

% Pointers to rows of x and y
L1=yr>=yrcgo & yr<=yrcsp; % calibration period
L2=isnan(y); % NaN data in y
L3=L1 & L2; % NaN data in calibration period in y

% Long term means and std devs of x and y
xmean=nanmean(x);
ymean=nanmean(y);
xstd=nanstd(x);
ystd=nanstd(y);

% Stdzd anomalies of x from its long term mean
x1= (x-xmean)/xstd;

% Synthetic y based on stdzd anomalies of x
y1 = ymean + x1*ystd;

% Replace calibration period NaNs in yrev with synthetic y
yrev(L3)=y1(L3);

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

function [ikey,isurr,rsurr,nsurr,errflg]=highr2(R,j1,j2,j3,mlap,rmin)
% highr2: highest correlation between tree-ring series; surrogate series
% CALL: [ikey,isurr,rsurr,nsurr,errflg]=highr2(R,j1,j2,j3,mlap,rmin);
%
% Meko 4-19-97
%
% Subfunction written for use with surrgt?.m functions that
% identify surrogate tree-ring series
%
%******************* IN *********
%
% R (mR x 4)r correlation, sample size, series 1 index, series 2 index
%    as obtained from r4sov.m
% j1 (1 x nsers)L  pointer to completed series
% j2 (1 x nsers)L  pointer to yet to be completed series
% j3 (1 x nsers)L  pointer to special (specified series)
% mlap (1 x 1)i  minimum number of acceptable overlap years for
%    a correlation coefficient
% rmin (1 x 1)r minimum acceptable correlation
%
%j
%****************** OUT *********************
%
% ikey (1 x 1)i  index of newly completed series
% isurr (1 x 1)i index of its surrogate
% rsurr (1 x 1)r  correlation coefficient between key and surrogate
% nsurr (1 x 1)i  sample size (# obs) for rsurr
% errflg (1 x 1)i error flag
%    0 = no problem
%    1 = no pairs of j1,j2 series in cols 3,4 of R, after
%      ruling out j3 series
%    2 suitable pairs exist, but none with at least mlap observations
%      overlap for the computed r
%    3 suitable pairs exist, and overlap OK for at least 1, but
%      none have high enough correl coef, as specified by rmin
%
%
%***************** NOTES **********************************
%
% The correlations and associated info, R, are assumed to 
% have been generated by r4sov.m.  R is therefore assumed not
% to have any NaN elements.  highr2.m checks anyway and reports an
% error if R contains any NaNs
%
% Typically, the j1,j2,j3 logical series correspond to pointers
% to time series:
% j1 -- 'completed' series.  Those either covering whole 
%   calibration period, or if not, those previously filled in
%   for the calibration period 
% j2 -- 'waiting' series.  Those not yet with surrogate data
%    generated for the calibration period
% j3 -- 'specified' series.  Series that the user does not want
%   correlations to determine surrogate for.  The user instead
%   has specified the surrogates based on other knowledge.

errflg=[]; %initialize error flag
rsurr=[];
nsurr=[];
isurr=[];
ikey=[];

% Check inputs
[mR,nR]=size(R);
if nR~=4,
   error('R must be 4-col matrix');
end
if any(any(isnan(R)));
   disp('R, as generated by r4sov.m, should not contain any NaNs')
   error('NaN elements in R');
end



Ltemp=[size(j1,1)==1 size(j2,1)==1 size(j3,1)==1];
if ~all(Ltemp);
   error('j1, j2, j3 must all be row vectors');
end
ntemp=size(j1,2);
Ltemp=[ntemp==size(j2,2) ntemp==size(j3,2)];
if ~all(Ltemp);
   error('j1,j2,j3 not all same col size');
end
Ltemp=[islogical(j1) islogical(j2) islogical(j3)];
if ~all(Ltemp);
   error('j1,j2,j3 not all logical');
end

% Check that j1, j2 mutually exclusive
Ltemp=j1 & j2;
if any(Ltemp),
   error('Some element is turned on in both j1 and j2');
end

% Check that have at least one j2 series and at least one
% j1 series, and that those series are not j3 series
Ltemp=j1 & ~j3;
if ~any(Ltemp);
   error('No j1 series, or none that are not also j3');
end
Ltemp=j2 & ~j3;
if ~any(Ltemp);
   error('No j2 series or none that are not also j3');
end



%**************************************************************
% Identify rows of R corresponding to pairs of series such
% that (1) one series is a j1 series but not a j3 series, and
% (2) the other series is a j2 series but not a j3 series, (3)
% the sample size for correlation is at least mlap observations,
% and (4) the correlation is at least rmin

% Acceptable sample size?
L1=R(:,2)>=mlap;

% Accepable minimum r
L2=R(:,1)>=rmin;

% Store cols 3 and 4 of R in a and b
a=R(:,3); % cv of indices, first of pair
b=R(:,4); % cv of indices, second of pair


% First or second series (col 3 or 4 of R) not a j3 series
if ~any(j3);
   L3a=logical(zeros(mR,1));
   L3b=L3a;
else
   numj3=sum(j3);
   j3f=find(j3); % a rv of indices of j3 series
   j3f=repmat(j3f,mR,1); % expand to a matrix, same row-size as R
   
   % Handle col 3 of R
   A=repmat(a,1,numj3); % expand to same size as j3f
   Ltemp=A==j3f;
   if numj3>1;
      L3a=(any(Ltemp'))';
   else
      L3a=Ltemp; % special case for vector instead of matrix
   end
   
   
   % Handle col 4 of R
   B=repmat(b,1,numj3); % expand to same size as j3f
   Ltemp=B==j3f;
   if numj3>1
      L3b=(any(Ltemp'))';
   else
      L3b=Ltemp;
   end
end


% Col 3 of R is a j1 series
numj1=sum(j1);
j1f=find(j1);
j1f=repmat(j1f,mR,1);

A=repmat(a,1,numj1);
Ltemp=A==j1f;
if numj1>1;
   L4a=(any(Ltemp'))';
else
   L4a=Ltemp;
end


% Col 4 of R is a j1 series
B=repmat(b,1,numj1);
Ltemp=B==j1f;
if numj1>1;
   L4b=(any(Ltemp'))';
else
   L4b=Ltemp;
end


% Col 3 of R is a j2 series
numj2=sum(j2);
j2f=find(j2);
j2f=repmat(j2f,mR,1);

A=repmat(a,1,numj2);
Ltemp=A==j2f;
if numj2>1;
   L5a=(any(Ltemp'))';
else
   L5a=Ltemp;
end


% Col 4 of R is a j2 series
B=repmat(b,1,numj2);
Ltemp=B==j2f;
if numj2>1;
   L5b=(any(Ltemp'))';
else
   L5b=Ltemp;
end


% Now we have these pointers to rows of R:
% L1 -- at least mlap observations for computation of R
% L2 -- r at least as large as rmin
% L3a -- j3 series in col 3
% L3b -- j3 series in col 4
% L4a  -- j1 series in col 3
% L4b -- j1 series in col 4
% L5a -- j2 series in col 3
% L5b -- j2 series in col 4


% One a j1 series, one a j2 series, and neither a j3 series
L6=((L4a | L4b) & (L5a | L5b)) & ~L3a & ~L3b;
if sum(L6)==0;
   disp('No j1,j2 pairs in R after excluding j3 series');
   pause(2);
   errflg=1;
   return
end

% Constraint L6 and: minimum overlap of mlap observations
L7=L6 & L1;
if sum(L7)==0;
   disp('Acceptable j1,j2 pairs, but not enought overlap in yr');
   pause(2);
   errflg=2;
   return
end

% Constraint L7 and: correlation at least rmin
L8=L7 & L2;
if sum(L2)==0;
   disp('Bombed out because correlation coef lower than rmin');
   pause(2)
   errflg=3
   return
end


%***************************************************
%
% L8 is the pointer to rows of R for which to find the max correl

f8=find(L8);  % row index in R of candidate rows

r1=R(L8,1); % cv of correlation coefficients
[rmax, imax]=max(r1);  % highest correlation coef, and row index to r1

irow = f8(imax); % row of R containing highest correlation
rsurr=R(irow,1);
nsurr=R(irow,2);
var3=R(irow,3);
var4=R(irow,4);

% Know that var3 and var4 are a j1,j2 pair, but do not know
% which is the j1 and which is the j2
if any(find(j1)==var3);
   ikey=var4;
   isurr=var3;
else
   ikey=var3;
   isurr=var4;
end
errflg=0;


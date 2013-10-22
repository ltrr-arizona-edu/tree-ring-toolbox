function [Y,J,D,ord,nump]=surrgt4 (X,YRS,JS,mlap,rmin,kopts)
% surrgt2.m: surrogate tree-ring indices  
% CALL: [Y,J,D,ord,nump]=surrgt4 (X,YRS,JS,mlap,rmin,kopts)
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
% Unlike surrgt2.m, surrgt4.m establishes heirarchy of daisy chain
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
% surrgt3a.m -- returns revised series given original series and surrogate series
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

         
         
         

         
   
   



x=x;





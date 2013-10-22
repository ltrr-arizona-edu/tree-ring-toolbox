function [W,N,J]=wgtdist1(E,L,dcrit,dsrch,nmax,kopts)
% [W,N,J]=wgtdist1(E,L,dcrit,dsrch,nmax,kopts)
%
% Inverse distance or inv-dist squared weights of gridpoints or tree sites on 
% stations for time series with time varying station makeup
%
% Meko 1-29-97
%
%******************** IN ARGS *****************************
%
% E (mE x nE)r  distances from mE gridpoints to nE stations. Get 
%	this previously using gcdist.m on the matrices of long, lat
% L (mL x nL)i  logical matrix telling for each of mL years
%		which of the nL stations have data (1) or missing
%		data (0).  nL==nE.
% dcrit(1 x 1)r threshold distance (km).  Any element of
%		E less than dcrit will be set to dcrit to avoid overemphasis
%		on very close stations
% dsrch (1 x 1)r search radius (km) around gridpoints. 
%		When kopts(1)==2, only the nearest up-to-nmax stations
%		within the search radius are used for the gridpoint. 
%		If kopts(1)==1, the search radius means nothing, as
%		the nearest nmax stations are used, regardless of 
%		whether in the search radius
% nmax (1 x 1)i maximum number of stations to use in weighting
%		the gridpoint.  Whether these need to be in the search
%		radius depends on kopts(1)
% kopts (1 x 2)i options:
%		kopts(1) ==1 nearest nmax stns, desregarding search radius
%				==2 nearest up-to-nmax stns within search radius
%		kopts(2) ==1 weighting method inverse distance
%		         ==2 ... inverse distance squared
%
%
%*******************  OUT ARGS *********************************
%
% W (mW x nW)r weights on stations for mW years.  nW==(mE*nE).
%		First nE cols of W apply to gridpoint 1, sencond nE cols
%		to gridpoint 2, etc.  Weighting method specified by
%		kopts(2)
% N (mN x nN)i  number of "predictor" stations for mN years
%		and nN gridpoints. Note that mN==mL and nN==mE.
% J (mJ x 2)i  pointer to cols of W telling start and end
%		columns applicable to each of mJ gridpoints.  Note 
%		that mJ == mE
%
%*** NOTES ******************************************************
%
% If kopts(1)==2 and no stations within search radius in any year, output
% arguments will be W all NaN and  N all zero
% 
% Quality control of input

[mE,nE]=size(E); % mE is number of gridpoints, nE number of stns
[mL,nL]=size(L); % mL is num of years, nL is number of stns
[m1,n1]=size(kopts);
if nE~=nL,
	error('E and L must have same col size');
end
if m1~=1 | n1~=2,
	error('kopts should be 1 x 2');
end

nyrs = mL; % number of years, or number of observations


% Size output matrices

a=NaN;
nW=mE*nE; % W will have this many columns
W=a(ones(nyrs,1),ones(nW,1)); % Weights
N=repmat(NaN,nyrs,1); % number of predictor stations per gridpoint

J=a(ones(mE,1),ones(2,1));

% Compute start, end col in W for each gridpoints weights
j1=(1:mE)';
jgo=1 + nE*(j1-1);
jsp=jgo+nE-1;
J=[jgo,jsp];

% If kopts(1)==2, substitute NaNs for any element of the
% distance matrix greater than the search radius
if kopts(1)==2;
	L1 = E>dsrch;
	sum1=sum(sum(L1));
	if sum1>0;
		E(L1)=a(ones(sum1,1),:);
	end
end

% Substitute dcrit for any elements of E less than dcrit
Ltemp=E<dcrit;
sumtemp=sum(sum(Ltemp));
if sumtemp>0;
	E(Ltemp)=dcrit(ones(sumtemp,1),:);
end


% Loop over the gridpoints
for n=1:mE;
	jgo=J(n,1);
	jsp=J(n,2);
   e=E(n,:); % rv of revised distances for this gridpoint -- with NaN for any station
   %  not within search radius (depending on kopts), and dcrit for any station
   % in the critical distance
   F = repmat(e,nyrs,1);    % F=e(ones(nyrs,1),:); % expand rv e into a matrix same 
   % size as L
   % Change elements of F corresponding to no data that year to NaN
	LN=~L; % 1 if no data for a station/year/month
   sum1=sum(sum(LN));  % Total number of NaNs in ts matrix of years by stations
   
   % If any NaNs in the years x stations matrix, change the corresponding values in 
   % the distance matrix to NaN
   if sum1>0;
      F(LN)= NaN;
   end
   
	% Error if kopts(1)==1 and any years have no acceptable stations
	L1=isnan(F); % 1 means a NaN in the years x stations distance matrix
	L2=all(L1'); % rv with 1 if a year has NaN at every station 
	if any(L2); % if any year has NaN at every station
      itemp=find(L2);
      % If disregarding search radius, should have stations "pulled in" for
      % every year.  Error if not.
      if kopts(1)==1; % if disregarding search radius
         disp(itemp);
         error('The above observations have no stations with data');
      end;
      
   end
   
   %**********************************************************
   % F is a matrix of distance for this gridpoint.  The rows of F are years.
   % the cols of F are stations.  This F is the distance from the gridpoint to
   % stations for all years.  Previous steps have substituted NaNs or dcrit as
   % needed in elements of F in accord with the settings of dcrit and dsrch if necessary. 
   % But want only the nearest nmax or fewer stations, depending on kopts.  
   % So want to rank elements in rows of F, then turn all except the lowest nmax
   % values to NaN
   
   Forig=F;
      
	[Fsort i1]=sort(F');
	[m2,n2]=size(Fsort); % m2 is number of station, n2 is number of years
   
   % Fsort now has dimensions stations x years,and is sorted from smallest to
   % largest distance. Want to cull the nearest nmax stations, corresponding to the
   % first nmax rows of Fsort.  Set the other rows to NaN,
   %%num1=n2*(m2-nmax); % this many elements of Fsort must be NaN'd
   %%Fsort((nmax+1):m2,:)=a(ones(m2-nmax,1),ones(n2,1));
   Fsort((nmax+1):m2,:)=NaN;  % I think this does the same thing as the prev 2 stmts
   
   % Fsort now has dimensions stations x years
	% Need to put elements of Fsort (now revised) into their original slots
	% in F.  To do this, convert matrices to column vectors, and convert
	% the row,col reference matrix i1 to col vector reference.  First build
	% cv of row indices
	j1=(1:n2)-1;
	j2=j1*m2; % a rv 
	J2=j2(ones(m2,1),:); % duped to a matrix
	i2=i1+J2;
	i3=i2(:);

	% Next string out Fsort and F
   Fsort=Fsort(:);
   F=F'; %  REVISION  ; % F is now stations by years
   F=F(:);
   % F and Fsort are now column vectors, strung out so that first rows are 
   % stations for first year, next rows are all stations for second year, etc
   
   % Recall that m2 is number of stations and n2 is number of years
	% Now substitute from the sorted and altered matrix, and reshape
	F(i3)=Fsort;
   F=reshape(F,m2,n2);
   
     
   % F is a matrix with as many rows as stations and as
   % many columns as years. Compute weights from the distances for the 
   % current gridpoint.
   if kopts(2)==1; % inverse distance
      G=1 ./ F;
   elseif kopts(2)==2;  % inverse distance squared
      G=1 ./ (F .* F);
   end;
   
   % Sum over rows (stations) to get the sum of weights for each year
   [sumw,nvals]=sumnan(G);
   
   % sumw is a row vector with length equal to the number of years. Dupe 
   % rows to to make SUMW, a 2-d matrix with number of rows equal to number
   % of stations.
   SUMW=repmat(sumw,m2,1);
   % OLD sumw(ones(m2,1),:);
   
   % Scale weights so that for each year they sum to 1.0 and 
   % Store scaled weights in multi-gridpoint matrices.  Use L3 and L4 to 
   % avoid 'divide by zero' warning when NaN in G divided by zero in SUMW
   L3=SUMW==0;
   L4=~isnan(G);
   L5=L3 & L4; 
   if any(any(L5)); 
      error('An element of G is non-NaN and the corresponding element of SUMW is zero');
   end;
   SUMW(L3)=NaN;
   
   % Compute weights and store in weighting matrix
   W(:,jgo:jsp)=(G ./ SUMW)';  % transpose makes this nyears by nstations
   N(:,n)=nvals';
end

function [W,N,J]=wgtdist2(E,L,dcrit,dsrch,nmax,kopts)
% [W,N,J]=wgtdist2(E,L,dcrit,dsrch,nmax,kopts)
%
% Inverse distance weights of stations on a target point.  
% Written for use with regcli3b.m in generating monthly climate
% series for a tree location or gridpoint.  Unlike wgtdist1.m, 
% wgtdist2.m is designed to work on only a single target point
%
% Meko 11-24-97
%
%******************** IN ARGS *****************************
%
% E (1 x nE)r  distances from the point to nE stations. 
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
%			(only option so far)
%
%
%*******************  OUT ARGS *********************************
%
% W (mW x nW)r weights on stations for mW years.  nW==nE.
%		Column j of W holds the weight on station j
%		Weighting method specified by kopts(2)
% N (mN x 1)i  number of "predictor" stations for mN years
% J (1 x 2)i  pointer to cols of W telling start and end
%		columns of weights
%
%*********************************************************

% Quality control of input

[mE,nE]=size(E); % mE is number of gridpoints, nE number of stns
[mL,nL]=size(L); % mL is num of years, nL is number of stns
[m1,n1]=size(kopts);
if nE~=nL,
	error('E and L must have same col size');
end
if mE ~=1;
	error('E must have row-size 1');
end
if m1~=1 | n1~=2,
	error('kopts should be 1 x 2');
end

nyrs = mL; % number of years, or number of observations


%-------------------- Size and allocate output matrices

a=NaN;
nW=nE; % W will have this many columns
W=a(ones(nyrs,1),ones(nW,1)); % Weights
N=a(ones(nyrs,1),:); % number of predictor stations per gridpoint
J=a(:,ones(2,1));


% ----------------Assign start, end col in W for weights
J=[1 nE];

% ---------------- Special treatment if require stations to be in dsrch
% If kopts(1)==2, user requires that station be within dsearch km
% of point to be considered a candidate for weighting.  Assign  NaNs
% to any element of the distance vector greater than the search radius
if kopts(1)==2;
	L1 = E>dsrch;
	sum1=sum(sum(L1));
	if sum1>0;
		E(L1)=NaN;
	end
end



%-------------- Substitute dcrit for any elements of E less than dcrit
Ltemp=E<dcrit;
sumtemp=sum(sum(Ltemp));
if sumtemp>0;
	E(Ltemp)=dcrit;
end


%----- Expand distance rv to matrix, and put NaNs where data missing
jgo=J(1,1);
jsp=J(1,2);
e= E; % rename distances for this gridpoint
% expand rv e into a matrix same size as L
F = e(ones(nyrs,1),:);
% Change elements of F corresponding to no data that year to NaN
LN=~L;
sum1=sum(sum(LN));
if sum1>0;
	F(LN)=NaN;
end

% Build cv logical flag to any years for which none of the nE stations have data
L1=isnan(F);
L2=(all(L1'));



%**********************************************************
% F contains distances from point to each station (col).
% Allowance has been made for dcrit and dsrch if necessary.  But
% want only the nearest nmax or fewer stations.  So want to rank
% elements in rows of F, then turn all except the lowest nmax values to
% NaN

[Fsort i1]=sort(F');
[m2,n2]=size(Fsort); % m2 stations, n2 years

% Set all except rows 1-to-nmax of Fsort to NaN
num1=n2*(m2-nmax); % this many elements of Fsort must be NaN'd
Fsort((nmax+1):m2,:)=NaN;

% Need to put elements of Fsort (now revised) into their original slots
% in F.  To do this, convert matrices to column vectors, and convert
% the row,col reference matrix i1 to col vector reference.  First build
% cv of row indices
j1=(1:n2)-1; % a rv of 0:(n2-1) , where n2 is number of years
j2=j1*m2; % a rv 
J2=j2(ones(m2,1),:); % duped to a matrix
i2=i1+J2;
i3=i2(:);

% Next string out Fsort and make strung-out vector with same number of 
% elements as F
Fsort=Fsort(:);
F=repmat(a,m2*n2,1);

% Now substitute from the sorted and altered matrix, and reshape
F(i3)=Fsort;
F=reshape(F,m2,n2);

% F now is m2 stations by n2 years, where all except the relevant
% distances (controlled by search criteria) are NaN
% Now compute the point-by-point inverse
% and sum over rows (stations) to get the weights and sum of weights for each year
G=1 ./ F;
[sumw,nvals]=sumnan(G); % nvals is number of valid stations in each year

% Set any value of sumw that is zero to NaN.  This because will use sumw
% elements as denominator in normalizing weights to sum to 1, and would
% get divide by zero error.
Lzero = sumw==0;
if any(Lzero);
   sumw(Lzero)=NaN;
end


% Make sum of weights 1.0
SUMW=sumw(ones(m2,1),:); % dupe the rv of sum of weights to same row size as G

% Compute and store weights
W(:,jgo:jsp)=(G ./ SUMW)';
N(:)=nvals';





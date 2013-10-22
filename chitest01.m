function [c,N,T,alpha]=chitest01(x,kopt,c1,ncrit);
% chitest01: chi square for normal distribution
% [c,N,T,alpha]=chitest01(x,kopt,c1,ncrit);
% Last revised  8-4-00
%
% Chi square test of data for normal distribution; parameters estimated from the data  
%
%*** INPUT ARGUMENTS
%
% x (m x 1)r data, length m
% c1 (1 x 1)i number of classes, or [], depending on kopt(1)
% kopt (1 x 2)i 
%	kopt(1)==1 automatic selection of number of classes (c1==[]) (see notes)
%				2 use c1 as number of classes (see notes)
%  kopt(2)==1 short-output mode
%				2 long-output mode
% ncrit(1 x 1)i  smallest allowable number of obs in any class
%
%*** OUTPUT ARGUMENTS
%
% c (1 x 1)i final number of classes used in the test 
% N (2 x c)i number of observed and expected in each class
% T (1 x 1)r test statistic
% alpha (1 x 1)r alpha level for significance of T
%
%
%*** REFERENCES
%
% Conover 
%
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED - none
%
%*** NOTES  ****************************************
%
% Number of classes.  If 
% Sheppard's correction
% Mean and std estimated from grouped data

mx=size(x,1); % length of x

% Compute max and min of data
xmax=max(x); xmin=min(x);

% Stretch data 1% on both ends
xrng = range(x);
xadd=xrng/100;
x1=xmin-xadd; x2=xmax+xadd;

% Initialize number of classes
switch kopt(1);
case 1;
   nc=2;
case 2;
   nc=c1-1;
end;

cbound=[];
n2=[];
kwh1 = 1; % while loop over number of classes
while kwh1 ==1; 
   nc=nc+1;
   xdiff=(x2-x1)/nc;
   cold=cbound;
   n2old=n2;
   cbound=x1+ (0:(nc))*xdiff; % class boundries
   B=repmat(cbound,mx,1);
   X=repmat(x,1,(nc+1));
   L1=X<B;
   s1=sum(L1);
   n2=diff(s1); % rv of # obs in each class
   
   if kopt(1)==2; % if number of classes specified
      if any(n2<ncrit);
         error(['Some class with fewer than ' int2str(ncrit) ' observations']);
      else;
         kwh1=0;
      end;
   else; % automatic selection of numnber of classes
      if any(n2<ncrit) | nc>9;
         if nc<=3;
            error(['Too few in a class and only three classes']);
         else;
            nc=nc-1;
            cbound=cold;
            n2=n2old;
            kwh1=0;
         end;
      end;
   end;
end; % kwh1

% Check that sum of obs in all classes equals sample size
if sum(n2)~=mx;
   error('Sum of members in classes not equal to mx');
end;

% Compute sample mean and standard deviation from grouped data
bmid=cbound+xdiff/2; % midpoints of classes
bmid(length(bmid))=[];
xbar=sum((bmid .* n2))/mx; % mean of grouped data
xbar1 = mean(x); % mean of ungrouped data
xvar = (((sum(n2 .* (bmid .^2)))/mx)-(xbar .^2)); % variance of grouped data
sheppard=(xdiff .^2)/12;  % Sheppard's correction quantity --
% to be subtracted from grouped variance before taking square root
xsd=sqrt(xvar-sheppard); % standard deviation, from grouped data, using 
% Sheppard's correction
xsd1 = std(x);

% OBSERVED NUMBER
O=[0 n2 0]; % tails on
nc=nc+2;

% HYPOTHESIZED PROBABILITIES

zbound = (cbound-xbar) ./ xsd; % Compute the standardized  lower class boundaries
   % of internal classes
F=normcdf(zbound,0,1); % theoretical distribution function at lower class boundaries
F=[0 F 1.0];  % prepare to compute theoretical probabilities for classes, including tails
pstar=diff(F); % probabilites;  end values apply to tails

% EXPECTED NUMBER
E=pstar * mx;

%******** COMBINING TAIL CLASSES IF TOO FEW EXPECTED OBS IN TAILS

%--- First on low end

EE=E(1:floor((nc-1)/2));
OO=O(1:floor((nc-1)/2));
a=cumsum(EE);
i=find(a>=ncrit);
if isempty(i);
   error('Cant get ncrit expected in low tails even combining');
end;
i=min(i);
if i==1; % no need to combine
else;
   aE=sum(EE(1:i));
   aO=sum(OO(1:i));
   E(i)=aE;
   E(1:(i-1))=[];
   O(i)=aO;
   O(1:(i-1))=[];
   nc = nc - (i-1);
end;

%--- Now high end

E=fliplr(E);
O=fliplr(O);
EE=E(1:floor((nc-1)/2));
OO=O(1:floor((nc-1)/2));
a=cumsum(EE);
i=find(a>=ncrit);
if isempty(i);
   error('Cant get ncrit expected in low tails even combining');
end;
i=min(i);
if i==1; % no need to combine
else;
   aE=sum(EE(1:i));
   aO=sum(OO(1:i));
   E(i)=aE;
   E(1:(i-1))=[];
   O(i)=aO;
   O(1:(i-1))=[];
   nc = nc - (i-1);
end;
E=fliplr(E);
O=fliplr(O);


% TEST STATISTIC
term1=O-E; % observed minus expected in each class
diffsq = term1 .^2; 
T=sum(diffsq ./ E);

% DEGREES OF FREEDOM
df = nc - 2 -1 ; % 2 parameters estimated

% Chi squared value
p=chi2cdf(T,df);
alpha=1-p;

% Arrange output
c=nc;
N=[O;E];



            
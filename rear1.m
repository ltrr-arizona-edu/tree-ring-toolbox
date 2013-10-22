% rear1.m   --- rearrange a monthly ppt array

% July, 1991, by D. Meko

% Main input array has monthly ppt for several stations.  One series after
% another as 13 cols and variable number of rows.  Output array is a 
% multivariate time series array with stations as sets of 12 columns.

% Also checks for consistency in number of missing values in each 
% station-month

%*****************  PRELIMS  ************************************

% Have an ascii file of trl-format monthly ppt series, one
%   	series after another, with one-line header
% Run FIL999.F to convert this to a numeric-only array suitable
%  	for matlab, while also replacing blank values with " 9999"
% In matlab, key in a string array (F) holding 6-char station names
% In matlab, save F and the approp outputs from FIL999.F in X.MAT

%****************  PRELOADS  ************************************

% X.MAT ... containing:

% A ... input ppt array  [?x13]
% B ... beg, end year for each series [?x2], as obtained from trl
% F ... names of stations  [?x6 str]  (six-character names)
% C ... missing values for each station/month [?x13]
%		(Produced by FILL999.F)
%*******************  OUTPUT  ************************************

% Y.MAT ... containing

% a ... station sequence [?x1]
% B ... as above
% b ... year sequence for rows of D
% D ... revamped ppt array
% F ... as above
% H ... row range in D for each ppt station [?x2]
% S ... col range in D ...

%*********************  VARIABLE LIST *****************************

% a  [m1x1] sequence number of series (e.g., [1 2 3...m1])
% A  [?x13] input monthly ppt array, first col the year
% b  [1x2] start, stop year of D
% B  [?x2] beg, end year for each ppt series
% b  [?x1] years for target array D
% C  [?x13] num mssg values, each series, each month; col 1 the ser #
% code...number to be filled in for missing values (e.g., 99.99)
%     code is "hard-coded".
% D [m1xn2] target output ppt array
% d [m1x1]  number of years in each component series of A
% e [1x2] two key years for echo table
% f [1x2]  e converted to subscripts relative to mgo of D
% F [m1x6 str] names of each series -- 6-character
% g [m1x1] additional "missing" years in D above that in A for each
% h [1xn2] # mssg values in each column of D
% H [m1x2] beg,end row index in D of each component series from A
%       These rows corresp to the ranges specified by B.
% I [m1x2]  beg, end row index in A of each component station's series
% J [m1x2]  cols of D corresp to Jan and July, for each station
% 	Each row a station.   This used in echo check for 2 key years.
% k1 ... Y/N for rescaling A by multiplying by 0.01.  This to avoid 
%     scaling more than once after initial load of X.MAT.
% K [m1x4]  Echo-check values of D for Jan, July of key years in e.
% 	K is also displayed for prtscr with heading.
% L [m1x12] number of mssg values in each station/month of D
%   Each row a station.  Each column a month.  Rearranged from h.
% M [m2xn2]  logic, 1 if corresp elmt of D missing value, 0 otherwise
% m1,n1 ... in B.  Number of stations equals m1.  n1=2.
% m2,n2 ... in D.  m2 is number of years (rows) in the final array D.
% mgo ... first year of D
% mstop ... last year of D
% Q [m1x12] duped cols of g
% R [m1x12] difference (L-K) minus C(:,2:13).  R holds the difference
%    between the computed # of missing values in D, and the
%    corresponding # in input array C.  The two should match.  An
%    error message is given if they do not match.
% S [m1x2] cols in D corresp to Jan and Dec of each station. 
%	Each row a station.
% t [1x4 logic] tests for existence of input arrays


%*********************  LOADS AND PREALLOCATES *********************
%
% Check that required input exists as variables
%
k1=[];

t = [exist('A')  exist('B') exist('F')  exist('C')];
if t~=[1 1 1 1]
		disp(' YOU FORGOT TO PRELOAD X.MAT ')
		keyboard;    % break to keyboard to allow loading X.MAT
end

% Preallocate space for main output array, and fill with 99.99

[m1,n1]=size(B);
code=99.99   % missing values code
mgo=min(B(:,1));   % earliest year of earliest series
mstop=max(B(:,2));  % latest year of latest series

m2=mstop-mgo+1;  % number of rows in target array
n2 =  m1*12   % number of cols in the target array

D = code (ones(m2,1),ones(n2,1));
b=[mgo:mstop]';

%*****************  SCALE A  *************************************

k1=input('SCALE A AGAIN? Y/N [N]: ','s')
if isempty(k1)
	k1='N';
end

if k1=='Y'
	A(:,2:13)= A(:,2:13)*0.01;
end

%**************** FILL THE TARGET ARRAY   *************************

% Various blocks of data from A will go into specific blocks of D
% Compute the indices ranges for the source and target arrays

% Compute row range in A for each series

d = B(:,2)-B(:,1)+1;    % number of years in each series
I = [cumsum(d)-d+1  cumsum(d)];  % beg and end index in A of each

% Compute target row index in D for each series from A

H = [B(:,1)+1-mgo(ones(m1,1),:)   B(:,2)+1-mgo(ones(m1,1),:)];

% Compute target col ranges for stations in D

S=[(1:12:m1*12-11)'  (12:12:m1*12)'];

%  Fill the target array

for i = 1:m1;
	D(H(i,1):H(i,2),S(i,1):S(i,2)) = A (I(i,1):I(i,2),2:13);
end

%********************  ECHO CHECK OF JAN, JULY FOR 2 KEY YEARS

e=input('RV OF TWO KEY YEARS FOR ECHO TABLE: ');

f=e-mgo+1;
a=(1:m1)';

J=[(a-1)*12+1  (a-1)*12+7];  % 

K=[(D(f(1),J(:,1)))'   (D(f(1),J(:,2)))'  (D(f(2),J(:,1)))' ... 
(D(f(2),J(:,2)))'];

home
clc
disp('    ECHO OF INPUT DATA')
disp(' ')
disp(['      stn           ',int2str(e(1)),'                ',int2str(e(2))]);
disp('                JAN      JULY      JAN      JULY')

disp(' ');
disp([a K])

%**********************   SAVE THE RESULTS IN A .MAT FILE ************

k2=input('SAVE RESULTS IN Y.MAT? Y/N  [Y]: ','s')
if isempty(k2), k2='Y';
end

if (k2=='Y')
	save Y D H S F  a b B
	clc
	disp('YOU HAVE SAVED Y,D,H,S,F,B, a and b in Y.MAT')
else
	disp('YOU HAVE NOT SAVED RESULTS IN Y.MAT')
end

%*************   TEST FOR CONSISTENCY OF MISSING VALUES   ************


g = m2(ones(m1,1),:)-(B(:,2)-B(:,1)+1);
Q = g(:,ones(12,1));   
M = D >= 99.90;

h = sum(M);

L=zeros(m1,12);
L(:)=h;
L=L';


R = (L-Q)-C(:,2:13);
if any(any(R))
	disp('WARNING:  R AND C DONT MATCH;  INCONSISTENT MSSG YEARS')
else
	disp('BONUS: R AND C MATCH; MSSG VALUES CONSISTENT')
end

function climcon1(pthfln1,txtin,datin)
% climcon1:  convert Great Sand Dunes formatted monthly data to 13-col  matrices
% CALL: climcon1
%
% Meko 5-8-98
%
%***************** IN ARGS
%
% pthfln1 (1 x ?)s  path and filename of input data file
% txtin{}  -- string input variables
%  {1} (1 x ?)s variable name for matrix holding input data (e.g., 'Xall')
%  {2} {1 x ?}s  file prefixes for output files
%          (e.g., {'pcp','mintmp'})
%
% datin{}  -- numeric input
%
% {1} (2 x 2)i  month, year or first data row and last data row in 
%        input matrix (e.g., [1 1951; 12 1997])
%
% {2} (1 x ?)r  subtract this amount from input variables in converting
% {3} (1 x ?)r  then multiply by this factor
%      note: [0 0 0 0 0], and [1 1 1 1 1], say,  amount to no conversion  
%
%
%**************** OUT ARGS
%
%******************* NOTES
%
% Written to convert spreadsheet derived climate data from Andy Lopez into .mat files
% with 13 cols in my regular monthly climate data format


%------------ Check input

if nargin~=3;
   error('Need 3 input args');
end

if ~ischar(pthfln1);
   error('pthfln1 should be char');
end

if size(txtin,1)~=1 | size(txtin,2)~=2 ; error('txtin should be 1 x 2 cell'); end;

if size(datin,1)~=1 | size(datin,2)~=3; 
   error('datin should be 1 x 3 cell'); 
end;
addamt = datin{2}; % will adjust variables by adding this in conversion
multamt = datin{3}; % will then multiply by this factor
nvbl = size(addamt,2);  % number of variables
if nvbl ~= size(multamt,2);
   error('datin{2} and datin{3} should be same size');
end

   
%-------------- Load file with matrix of input data
eval(['load ' pthfln1]);
eval(['Xall = ' txtin{1} ';']);



%--------------- Check col size for consistency
[m1,n1]=size(Xall);
if n1 ~=nvbl;
   error('column size of input data matrix inconsistent with datin{2} and {3}');
end

%---------------- Augment rows of Xall if needed to begin in Jan and end in Dec
mnyr = datin{1};
if mnyr(1,1)~=1;
   monadd1 = mnyr(1,1)-1;
   Xall = [repmat(NaN,monadd1,nvbl);  Xall];
end
if mnyr(2,1)~=12;
   monadd2 = 12-mnyr(2,1);
   Xall = [Xall; repmat(NaN,monadd2,nvbl1)];
end
[m2,n2]=size(Xall);


%-------------- Check that number of rows of Xall equals 12  times number of years
nyr = mnyr(2,2)-mnyr(1,2)+1;
nOK = nyr * 12;
if m2 ~= nOK;
   error('Number of rows in (augmented) Xall does not equal 12 times compute n of yr');
end
disp('here');


%------------ Loop over variables
yr = (mnyr(1,2):mnyr(2,2))'; % year vector for X in output files

for n = 1:nvbl;
   x = Xall(:,n);
   X = (reshape(x,12,nyr))';
   X = (X+addamt(n)) * multamt(n); % convert units if necessary
   X = [yr X];
   flnout = txtin{2}(n);
   eval(['save ' flnout{1} ' X;']);
end

   
   
   
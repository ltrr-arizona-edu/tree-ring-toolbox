function Y=daycon2(X,kopt)
% daycon2:  Convert a BOR spreadsheet format daily flow series to a daily tsm
% Y=daycon2(X,kopt);
%
% Convert a BOR spreadsheet format daily flow series to a daily tsm
%
%*** INPUT
% X input matrix of daily flow data;  for nits and format see kopt(1) and notes
% kopt (1 x 1)i options
%    kopt(1) type of input file (see notes)
%       ==1 BOR Clear Creek format (originally read in with readwk1.m)
%
% Y (? x 4)r output matrix
%   col 1:  year (calendar)
%   col 2:  month (1-12)
%   col 3:  day  (1-31 for Jan, 1-29 for Feb, 1-31 for Mar, etc); Feb 29 data 
%           for non-leap years interpolated (see notes)
%
%*** OUTPUT
%
% Y (? x 5)r daily time series matrix.  Col are:
%   1-year (calendar)
%   2-month (jan=1)
%   3-day of month (use 29 for Feb)
%   4-day number (1 to 366 within calendar year)
%   5-daily ave flow (cfs)
%     (see notes)
%
%*** NOTES
%
% Linear interpolation of missing data.   If some days have missing data, flagged
% as indicated above, or are days of month not existing (e.g. Feb 31), flow data
% are linearly interpolated from the two flanking values.
%
% Units.  I convert whatever the input units are to cfs
%
% Output matrix.  Output matrix has 4 numeric columns:
%   1-year
%   2-month
%   3-day of month
%   4-sequential day (1-366)
%   5-ave daily values of flow (cfs)
% Each year has 366 rows, matching the number of days in a non-leap year. Feb
% 29 on non-leap years is set to NaN
%
% Subsequent use in n-day minimum flow calculation.  Functions I write to
% do this type of analyis will automatically handle NaN February 29's. 
%
% kopt(1)==type of input file.
% ==1 Clear Creek data was an xls spreadsheet that I used Excel to save as a .wk1
%     file.  Then used readwk1.m to read in and save in a matrix. Contents:
%     Col 1== water year
%     Col 2== month
%     Col 3== monthly total in acre-ft
%     Col 4== monthly ave in cfs (sum of daily cfs values)
%     Col 4-15 == daily ave flow in cfs
%        Conversion used by BOR in making col 3 from col 4 is  cfs = 1.9835 AF/day



%---  CHECK INPUT DATA

% Clear Creek format

if kopt(1)==1;
     
   % column size must be 35 and row size a multiple of 12
   [mX,nX]=size(X);
   if nX~=35;
      error('X must have 35 cols');
   end;
   if rem(mX,12) ~=0;
      error('Row size of X must be multiple of 12');
   else;
      nyr1 = mX/12; % number of water years
   end;
   
   %  First month must be Oct (10), last Sep (9)
   if X(1,2)~=10 | X(mX,2)~=9;
      error('First month of X must be 10 and last month 9');
   end;
   
   % Prepend 9 months of NaN data (Jan-Sep), and append 3 months (Oct-Dec)
   X1=[repmat(X(1,1),9,1)  (1:9)'  repmat(NaN,9,33)];
   X2=[repmat(X(mX,1),3,1)  (10:12)'  repmat(NaN,3,33)];
   X=[X1; X; X2];
   mX=size(X,1);
   if mX ~= nyr1*12 +12;
      error('Row size of X wrong after prepending and appending data');
   end;
   nyrX=mX/12; clear nyr1 X1 X2 
   
   
   yrgo = X(1,1)-1;
   yrsp = yrgo+nyrX-1;
   
   % Build calendar year vector
   i1 = 0:1:(nyrX-1);
   I1=repmat(i1,366,1);
   J1= yrgo + I1;
   yr = J1(:); % year vector
   clear i1 I1 J1 ;
   
   % Build month vector;
   j=[repmat(1,31,1); repmat(2,29,1);  repmat(3,31,1);  repmat(4,30,1); ...
         repmat(5,31,1);  repmat(6,30,1);  repmat(7,31,1);  repmat(8,31,1);...
         repmat(9,30,1);  repmat(10,31,1);  repmat(11,30,1); repmat(12,31,1)];
   J=repmat(j,1,nyrX);
   mon = J(:);
   clear J;
   
   % Build day of month vector
   
   j=[1:31   1:29   1:31 1:30 1:31 1:30 1:31 1:31 1:30 1:31 1:30 1:31]';
   J=repmat(j,1,nyrX);
   daymon= J(:);
   clear J;

   
   % Build day vector
   j=(1:366)';
   J=repmat(j,1,nyrX);
   day = J(:);
     
   
   % Allocate Y, and store year, month, day-of-month, julian day
   mY = nyrX*366; 
   Y = repmat(NaN,mY,5);
   Y(:,1:4)=[yr mon daymon day];
   
      
   % --------- Cull the monthly data and put in slots in Y
   
   Z = X(:,5:35); % store daily-data-only part of X as Z
   
   % String out daily data as vector
   Z=Z';
   z=Z(:); 
   
   % Build row-index vector to elements of z to be deleted. These would be:
   n31=31*12; % number of days in year in which each month has 31 days
   i1 = (1*31)+[30 31]; % Feb 30,31
   i2 = (3*31)+[31] % Apr 31
   i3 = (5*31)+31; % june 31
   i4 = (8*31)+31; % sep 31
   i5 = (10*31)+31; % Nov 31
   iout = [i1 i2 i3 i4 i5]';
   nout=length(iout);
   clear i1 i2 i3 i4 i5;
   
   i1=[0: n31:  (nyrX-1)*n31];
   I1 = repmat(i1,nout,1);;
   IOUT = repmat(iout,1,nyrX);
   I3 = IOUT+I1;
   I3=I3(:);
   clear I1 IOUT i1;
   
   % Remove "bogus" days, leaving 366 per year
   z(I3)=[];
       
   % Store daily data 
   Y(:,5)=z;
   clear z I3
   
   % Put NaN in non-leap year Feb 29
   L1 = Y(:,2)==2; % feb
   L2 = Y(:,3)==29; % 29th
   L3 =  rem(Y(:,1),4) ~= 0;
   L4 = L1 & L2 & L3;
   Y(L4,5)=NaN;
      
      
   else;
end;
disp('here');
   
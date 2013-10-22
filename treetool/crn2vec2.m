function [x,s,yr]=crn2vec2(pf1)
% crn2vec2:  .crn file to  column vectors of index, sample size, and year
% [x,s,yr]=crn2vec2(pf1);
% Last revised 6-29-99
%
% Reads a ".crn" file of tree-ring indices in ITRDB format and converts
% the indices and associated data into vectors that can be used in Matlab.
% Input file must be formatted as in ITRDB requirements for crn2vec2 to work.
% 
%*** IN *********************************************
%
% pf1 (1 x ?)s  path and file name of .crn file (e.g., 'c:\work\cityrks.crn')
%
%
%*** OUT ****************************************
%
% x (mx x 1)r   tree-ring index, mx years
% s (mx x 1)r   sample size (number of cores) in each year
% yr (mx x 1)i  year vector for x and s
%
%*** REFERENCES -- None
%
%*** UW FUNCTIONS CALLED -- None
%
%*** TOOLBOXES NEEDED -- None
%
%*** NOTES ****************************
%
% A related function, crn2vec1.m (not in tree-ring toolbox) handles many chronologies
% at a time

a=NaN; % will use this later
xa = a(ones(10,1),:);  % initialize 10 element NaN vector

% Open .crn file
fid1=fopen(pf1,'r');

% Skip first header line
c = fgetl(fid1);

% Get year info off second header line
c = fgetl(fid1);
yrgo=str2double(c(67:71));
yrsp=str2double(c(72:76));
yr = (yrgo:yrsp)'; % year vector for output
nyrs = yrsp-yrgo+1; % number of years of data

% Compute expected number of "decade lines" of data.  Note special case handling
% of negative years.  These occur when series are BC  and the file does not
% use the "add 8000" convention
ii1=yrgo-rem(yrgo,10); % first 'decade'
% Handle the 'negative' starting year case (a la El Malpais)
if rem(yrgo,10)<0;
   ii1=ii1-10; % start year in, say, -136, gives first decade as -140, not -130
end
ii2 = yrsp - rem(yrsp,10);  % ending decade
% Handle 'negative' end year
if rem(yrsp,10)<0;
   ii2=ii2-10;
end
nlines = round(ii2/10)-round(ii1/10) + 1;  % number of 'decade lines' of data

% Compute number of values to skip in reading first decade-line.
% Will differ depending on whether starting year is positive or negative.
% Negative first year example is El Malpais, NM
if yrgo>=0;
   nskip=rem(yrgo,10);
   % Let's say data starts in 1652. Code would say skip first 2 values on "1650" line
else; % negative first year
   nskip = 10 - abs(rem(yrgo,10)); 
   if nskip==10; nskip=0; end;
   % if yrgo is -139, this says skip one value on the "-140" decade line
end

% Compute number of values to read on last decade line
% Again, different treatment if yrsp is negative
if yrsp>=0;
   nlast = rem(yrsp,10)+1;  
   % For example, if yrsp is 1686, says read 7 values -- 1680 to 1686
else;  % yrsp is negative
   nlast = 10 - abs(rem(yrsp,10)) + 1;
   % So if yrsp is -59, says read 2 values (-60 and -59) off last line
   % Or if yrsp is -51, says read 10 values (-60,-59,...,-51) of last line
end


% Compute number of expected years of data
nyrs = yrsp-yrgo+1;

% Skip third header line
c=fgetl(fid1);

x=a(ones(nyrs,1),:);
s=x;

% Compute row slots for each within x,s for each line of data
irow = a(ones(nlines,1),ones(2,1));
irow(:,1) = 1 + 10 * ((0:(nlines-1))');
irow(2:nlines,1) = irow(2:nlines,1) - nskip;
irow(:,2) = irow(:,1)+9;
irow(1,2)=irow(1,2)-nskip;
irow(nlines,2) = nyrs;

% Initialize cols for data and sample size form most lines
igo1 = [11:7:74];  isp1 = igo1 + 3;  % start, end cols for data
igo2= isp1+1;  isp2 = igo2+2;  % start, end cols for sample size


for n= 1:nlines
   irgo =irow(n,1);
   irsp= irow(n,2);
   c=fgetl(fid1);
   % Compute start and end positions for the data values and sample size
   ion1=igo1;  ioff1=isp1; ion2=igo2;  ioff2=isp2; nvals=10;
   x1 = xa; s1=xa;  nvals=10;
   if n==1;
      ion1(1:nskip)=[]; ioff1(1:nskip)=[]; ion2(1:nskip)=[]; ioff2(1:nskip)=[];
      nvals=length(ion1);
      
      
   end
   if n==nlines;
      ion1=ion1(1:nlast);  ioff1=ioff1(1:nlast); ion2=ion2(1:nlast); ioff2=ioff2(1:nlast);
      nvals=length(ion1);
   end
   
   x1=xa(1:nvals);
   s1=xa(1:nvals);
   for m = 1:nvals;
      cdat=c(ion1(m):ioff1(m));
      sdat=c(ion2(m):ioff2(m));
      x1(m)=str2double(cdat);
      s1(m)=str2double(sdat);
   end
   x(irgo:irsp)=x1;
   s(irgo:irsp)=s1;
      

   
end

% Subtract 8000 from years vector if ending year greater than 2100 
if max(yr)>2100;
	yr = yr - 8000;
end


fclose(fid1);

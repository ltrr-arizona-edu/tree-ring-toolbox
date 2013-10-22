function decade02(x,yrx,s,nmax,k)
% decade01:  ascii .txt file of decade-formatted time series, append version
% decade01(x,yrx,s);
% Last revised 6-24-01
%
% First use was to build an appendix listing reconstructed annual flow of Sacramento River
%
%*** INPUT
% 
% x (mx x 1)r  time series vector to be listed
% yrx (mx x 1)i year vector for x
% s{} cell of string input
%   {1} (1 x ?)s title (e.g., 'Appendix 1.  Reconstructed Flow')
%   {2} (1 x ?)s format for line with row and 10 values
%   {3} (1 x ?)s path\filename of output .txt file
%   {4} (1 x ?)s header line with 0 1 2 ...9
% nmax (1 x 1)i maximum number of data lines per page
% k(1 x 1)i flag for series
%  k==1 first
%  k==2 last
%  k==3 intermediate
%
%*** OUTPUT
%
% No args.
% A .txt file with file in decade format is produced
%
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Typical Use.  Run the function.  Open the output file in MS word.  Change formatting as needed.  Add Opening 
% header line :  Year  0    1     3   ...   9  and put at top of each page
% Add appendix footer
% Save .doc file

% Unload

 strtit=s{1};
 fmt1=s{2};
 pfout=s{3};
 shead=s{4};

% Strip off NaN years
L=isnan(x);
x(L)=[];
yrx(L)=[];
d = diff(yrx);
if~all(all(d==1));
    error('Series has an internal NaN');
end;


% Calculate number of leading years
yrgo1 = yrx(1);
if yrgo1>=0;
    nlead = rem(yrgo1,10);
else;
    nlead = rem(yrgo1,10);
    if nlead~=0;
        nlead=10-abs(nlead);
    end;
end;


% Calculate number of trailing years
yrsp1 = max(yrx);
if yrsp1>=0;
    ntrail = 10-rem((yrsp1+1),10);
else;
    ntrail = rem((yrsp1+1),10);
    if ntrail~=0;
        ntrail=abs(ntrail);
    end;
end;


% TACK ON LEADING AND TRAILING YEARS

if nlead~=0;
    yrgo2=yrgo1-nlead;
    x = [repmat(NaN,nlead,1); x];
else;
    yrgo2=yrgo1;
end;
if ntrail~=0;
    yrsp2=yrsp1+ntrail;
    x = [x ; repmat(NaN,ntrail,1)];
else;
    yrsp2=yrsp1;
end
yrx = (yrgo2:yrsp2)';

if length(yrx)~=length(x);
    error('x and yrx not same length');
end;

if rem(length(x),10)~=0;
    error('padded length of x not multiple of 10');
end;


% OUTPUT

if k==1;
    fid1=fopen(pfout,'w');
else;
    fid1=fopen(pfout,'a');
end;

% title and skipped line
if k>1;
    fprintf(fid1,'\f');
end;

    fprintf(fid1,'%s\n',strtit);
    fprintf(fid1,'%s\n',blanks(10));
    fprintf(fid1,'%s\n',shead);
    fprintf(fid1,'%s\n',blanks(10));
    
    
nlines = length(x)/10;
m=0;
for n = 1:nlines;
    m=m+1;
    jgo = n*10-9;
    jsp=jgo+9;
    year = yrx(jgo);
    y = (x(jgo:jsp));
    z = [year ; y];
    
    if m>nmax;
        m=0;
        fprintf(fid1,'\f');
        fprintf(fid1,'%s\n',strtit);
        fprintf(fid1,'%s\n',blanks(10));
        fprintf(fid1,'%s\n',shead);
        fprintf(fid1,'%s\n',blanks(10)); 
        fprintf(fid1,'%s\n',blanks(10));
    end;
    
        
    
    fprintf(fid1,fmt1,z);
end;

fclose(fid1);









function y=tsv2rwl(nm,x,yrx)
% tsv2rwl:  time series vector to rwl-segment
% y=tsv2rwl(nm,x,yrx);
% Last revised 2006-08-09
%
% Convert a time series of ring-width measurements to a character-block suitable for building a .rwl file.
%
%*** INPUT
%
% nm (1 x ?)s id of series -- max of 8 characters, and usually 6
% x (mx x 1)r  time series of measurements, assumed to be in hundredths of mm
% yrx (mx x 1)r  year vector for x
%
%
%*** NOTES
%
% cols 9-12 reserved for the year.  If years extend  back earlier than AD -999, cols 8-12 reserved
% cols 13-18, 19-24, ... 67-72  reserved for the ring width values for 10 years per line
% Rows other than first and last all have 10 values and begin with an even-decade year (e.g., 1710)
% Year at beginning of row is year of first value in that row
% First row -- data cols left justified
% Last row -- likewise, but with additional oddity of -9999 as the "data" value following last 
% valid measured width.  If series ends in a 9 year (e.g., 1999), last row has only the -9999 data
%
% Rev2006-08-09:  Changed the output rwl file to have a -9999 instead of a 999 as the dummy value indicating
% the previous year had the last valid measurement.  This seems consistent with the newer measuring programs.


% Check x
[mx,nx]=size(x);
if nx~=1;
    error('x must be col vector');
end;
if size(yrx,2)~=1 | size(yrx,1)~=mx;
    error('yrx must be same length as x');
end;


if intnan(x)==1;
    smess=['Warning: ' nm ' has internal NaNs'];
    uiwait(msgbox(smess,'Warning','modal'));
end;


% slap on trailing -9999
x=[x; -9999];
yrx=[yrx; max(yrx)+1];

% Check nm
if ~isstr(nm);
    error('nm must be string');
end;
if length(nm)>8 | (length(nm)==8  & yrx(1)<-999);
    error('nm must be fewer than 9 chars, and fewer than 8 if first year earlier than AD -999');
end;


% Compute number of full (10-value) lines
j10 = find(rem(yrx,10)==0);
if length(x)<10;
    nmid=0;
    yr10=[];
else;
    nmid = length(j10);
    y=[];
    yr10=yrx(j10); % years starting rows in middle block
end;


% Store leading years
if nmid==0; % series so short that no decade years
    nlead = length(x);
    xlead=x;
    yrgolead=yrx(1);
elseif yrx(1)<yr10(1);
    nlead = yr10(1)-yrx(1);
    xlead = x(1:nlead);
    yrgolead = yrx(1);
else;
    nlead=0;
    xlead=[];
    yrgolead=[];
end;

% Store trailing years
yrlast = max(yrx);
if nmid==0;
    ntrail=0;
else;
    ntrail = yrlast-max(yr10)+1;
    if ntrail==10; % no trailing valies
        ntrail=0;
        % yr10 unchanged
    else;
        yrgotrail = max(yr10);
        L=yrx>=yr10(end) & yrx<=yrlast;
        xtrail=x(L);
        ntrail=length(xtrail);
        yrfinish=yr10(end);
        yr10(end)=[];
    end;
end;

if isempty(yr10); % series so short that no middle block
    flagmid=1;
else;
    flagmid=0;
end;


if flagmid==0; % if a middle block exists
    
    % Middle block
    L = yrx>=yr10(1) & yrx<=(yr10(end)+9);
    xmid=x(L);
    
    % Reshape the middle block
    ksize=length(xmid);
    if rem(ksize,10)~=0;
        error('Not divisible by 10');
    end;
    nrow = ksize/10; % number of middle row
    xmid = (reshape(xmid,10,nrow))';
    
    
    
    % DATA BLOCK
    
    % Middle
    M=num2str(xmid,'%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f\n');
    
    % Change 3-19-04a:  for Matlab 7.0, appears that \n leaves a trailing blank col in M; get rid of if
       if all(isspace(M(:,end)));
        M(:,end)=[];
    end;
    % Change 3-19-04a - end
    
    [mM,nM]=size(M);
    npad = 60-nM;
    if npad>0;
        M = [repmat(blanks(npad),mM,1) M];
    end;
    
end; % if flagmid==0; % if a middle block exists


% Prepend leading line
if nlead>0;
    slead =num2str(xlead','%6.0f');
    nneed = (6*nlead)-length(slead);
    if nneed>=0;
        slead = [blanks(nneed) slead];
    end;
end;
if flagmid==0 & nlead>0;
    M=char(slead, M);
elseif flagmid==0 & nlead==0;
    % M=M;
else;
    M=slead;
end


% Append trailing line
if ntrail>0;
    strail =num2str(xtrail','%6.0f');
    nneed = (6*ntrail)-length(strail);
    if nneed>=0;
        strail = [blanks(nneed) strail];
    end;
    M=char(M,strail);
end;


   
% Years col

if nlead>0;
    if flagmid==0;
        yr10=[yrx(1) ; yr10];
    else;
        yr10=yrx(1);
    end;
end;
if ntrail>0;
    yr10=[yr10; yrfinish];
end;
if min(yrx)<-999;
    s = num2str(yr10,'%5.0f');
    jmax = 7;
else;
    s = num2str(yr10,'%4.0f');
    jmax=8;
end;
    

% ID column

nlines=size(s,1);
lenyear = size(s,2);
if lenyear==4;
    % Year will be in cols 9-12; reserve 1-8 for id
    lenid=8;
elseif lenyear==5;
    lenid=7;
elseif lenyear<4; % the char matrix of years needs to be padded on left to width 4
    lenid=8;
    nshort = 4-lenyear;
    s=[repmat(blanks(nshort),nlines,1)  s];
end;
clear nshort;
    

if length(nm)<lenid;
    nadd= lenid-length(nm);
    nm=[nm blanks(nadd)];
else;
    nm=nm(1:lenid);
end;
id = repmat(nm,nlines,1);

% COMBINE ID AND UEAR
B=[id s];

% COMBINE WITH DATA
y = [B M];

if flagmid==0 & size(y,2)~=72;
    error('y col size must be 72');
end;










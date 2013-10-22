function [x,person,when]=rwread5(pfinput,kdigits);
% rwread5:  read a .rw-format ring-width file
% [x,person,when]=rwread5(id,suff);
% Last revised 2006-08-09
%
% Read a .rw-format ring-width file
%
%
%*** INPUT
%
% pfinput (1 x ?)s   path/filename of input  <'c:\import\pdf02a.eww'>
% kdigits:  ==1  input rw data is hundredths of mm
%           ==2  input rw data is thousandths of mm
%
%*** OUTPUT
%
% x (mx x 2)r  year and ring-width measurements, in mm
% person (1 x 3)s or (1 x 2)s  initials of measurer
% when (1 x 6)i  date measurements completed or last meaurements done ( year month day hour minute second)
%
%*** REFEREBNCES -- NONE
%*** UW FUNCTIONS CALLED 
% deblankb
%
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES 
%
% Rev2006-08-09: to accept new argument kdigits and optionally consider input rw file as in units of 100ths or
% 1000ths of mm


%-- TRIM OFF TRAILING 999 or -9999


%--- VARIABLE NUMBER OF IN ARGS

if nargin<2;
    kdigits=1; % assume data in 100ths of mm
else
    L=kdigits==[1 2];
    if ~any(L);
        error('kdigits in rwread5.m must be 1 or 2');
    end;
end

%--- SCALING OF INPUT

switch kdigits;
    case 1;
        fscale=100;
    case 2;
        fscale=1000;
    otherwise
end



S= textread(pfinput,'%s','delimiter','\n','whitespace','');

% Check if last element 999 or -9999 and trim it off 
send = S{end};
send = deblankb(send); % deblank on left and right
if strcmp(send,'999') | strcmp(send,'-9999');
    S(end)=[];
else;
    S(end)=[] ; % trim off the eof line, which sometimes is after the 999 or -9999 line
    send=S{end};
    send = deblankb(send); % deblank on left and right
    if strcmp(send,'999') | strcmp(send,'-9999');
        S(end)=[];
    else;
        S(end)=[] ; % trim off the eof line
        send=S{end};
        send = deblankb(send); % deblank on left and righ
        if strcmp(send,'999') | strcmp(send,'-9999');
            S(end)=[];
        else;
            error(['No 999 or -9999 in last threee lines of ' pfinput]);
        end;
    end;
end;


%-- STORE INITIALS, DATE OF MEASURING, AND FIRST YEAR
s1=S{1};
s2=S{2};
s3=S{3};


%-- CHECK INITIALS -- s1

s1=deblankb(s1); % deblank flanks
if ~all(isletter(s1));
    error([s1 ' not all letters']);
else;
    person=s1;
    S(1)=[]; % Strip measurer line off 
end;


%-- CHECK DATE

% Shoul have two slashes
k=strfind(s2,'/');
if length(k)~=2;
    error(['Date line ' s2 ' does not have two slashes']);
end;

% strip trailing blanks, then leading blanks
s2=deblankb(s2);
s2=fliplr(deblank(fliplr(s2)));
klen = length(s2); % length of date

% Get month, day, year parts of date
k=strfind(s2,'/');
if length(k)~=2;
    error(['Date line ' s2 ' does not have two slashes']);
end;

% -- Check validity of fields 
mon = str2num(s2(1:(k(1)-1)));
if mon<1 |  mon>12;
    error (['month ' mon ' invalid']);
end;

day = str2num(s2((k(1)+1):(k(2)-1)));
if day<1 | day>31;
    error (['day ' day  ' invalid']);
end;

yr = str2num(s2((k(2)+1):klen));
[yrnow,dum1,dum2,dum3,dum4,dum5] = datevec(date); % yrnow is current year;
if yr<0 |  yr >yrnow;
    error(['Measurement Year ' str2num(yr) ' impossible -- negative, or in future']);
end;
if yr<30;
    yr=yr+2000;
    dformat=2;
elseif yr <100;
    yr=yr+1900;  % in case 4/23/80 type coding of year
    dformat=2;
else;
    yr=yr;
    dformat=23;
end;

%--- STORAGE FOR  DATE IN DATE FORMAT
d=[NaN NaN NaN  0 0 0];
d(1)=yr;
d(2)=mon;
d(3)=day;
dwhen=d;

S(1)=[]; % strip date of measuring line off



% Doublecheck computed date against original in m/d/yr format
sdate = datestr(d,dformat);
if ~strcmp(s2,sdate);
    error(['Dates ' sdate ' and ' s2 ' not same']);
end;
when=sdate; % return string date
when=datestr(datenum(when),23); % so that year in 4-digit form 

%--- MAKE YEAR VECTOR

yrgo = str2num(s3);
if isempty(yrgo);
    error(['Start year ' s3 ' invalid']);
end;
S(1)=[]; % strip off start year
nyr = size(S,1);
yrsp = yrgo+nyr-1;
yr = (yrgo:yrsp)';


%--- CONVERT MEASUREMENTS

x=str2num(char(S));
if isempty(x);
    error('Could not convert measurements to number');
end;
x = x/fscale;  % convert from 100ths or 1000ths of mm to mm


%--- COMBINE YEAR WITH MEASUREMENT

x=[yr x];











    
    

   





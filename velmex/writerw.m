function writerw(x,person,when,yrfirst,fn,suffx,path1)
% writerw: write rin width data as a .rw, .eww or .lww file
% writerw((x,person,when,yrfirst,fn,suffx);
% Last revised 2006-08-09
%
% Writes ring-width measurements as Bannister-machine-compatible .rw files 
%
%
%*** INPUT
%
% x (mx x 1)r ring-width measurements
% person (1 x ?)s   initials of measurerer (2 or 3 chars)
% when (1 x ?)s   string-format date of measurement (e.g.,  21-Sept'01)
% yrfirst (1 x 1)i  first year of measurements
% fn (1 x ?)s  filename of output .rw, .eww or .lww file
% suffx (1 x ?)s  suffix of output file (e.g., .rw, .eww or .lww
% path1(1 x ?)s director to to write file.  End path with a \
%
%
%*** OUTPUT
%
% No arguments. 
% A .rw, .eww or .lww file is written
%
%*** NOTES
%
% Units.  Units assumed to have been converted to hundredths of mm before call to rwc2rw
% Output format.  Output is ascii file with line:
%   1) initials
%   2) date
%   3) year of first measurement
%   4:?)  one value per line, hundredths of mm, left justified in column 2
%   999  end indicator
%
% Carriage return/linefeed.  Note that fprintf use \r\n to produce the car retn plus line feed
%
% Rev2006-08-09:  to write -9999 rather than 999 as end terminator
% Rev2011-04-04:  to handle path and file naming problems in linux

% Check x;
[mx,nx]=size(x);
if nx~=1 | mx<2;
    error('x must be col vector of length 2 or greater');
end
if ~isnumeric(x);
    error('x must be numberic');
end;
if any(x<0) | any(x>50000);
    error('x in writerw.m cannot be negative and cannot be greater than 50000 (50 mm or 500 mm, depending on units)');
end;


% Person
if ~isstr(person);
    error('person must be string');
end;
[mtemp,ntemp]=size(person);
if mtemp~=1;
    error('person must be row-string');
end;
if ntemp>3;
    error('Person must be 3 chars or fewer');
end;


% When
if ~isstr(when);
    error('when must be string');
end;
[mtemp,ntemp]=size(when);
if mtemp~=1;
    error('when must be row-string');
end;
if ntemp>11;
    error('When must be 11 chars or fewer');
end;
% Convert date to day/month/year -- e.g., want 9/21/2001 insteat of 21-Sept-01
dwhen=datestr(datenum(when),23);
yrcurr=datevec(dwhen);
yrcurr=yrcurr(1); % store the "year" of measurment for later check that not earlier than last measured year of data


% yrfirst
[mtemp,ntemp]=size(yrfirst);
if ~isnumeric(yrfirst) | (mtemp ~=1 & ntemp~=1);
    error('yrfirst must be 1 x 1 numeric scalar');
end;
% yrfirst and length of x consistent with current uear
yrlast = yrfirst+length(x)-1; % last year of measurements
if yrlast>yrcurr;
    error(['Last measured year (' int2str(yrlast) ') is more recent than date measurements entered or last edited']);
end;
 

% fn -- filename
if ~isstr(fn);
    error('fn not a string');
end;
if size(fn,1)~=1;
    error('fn has row dimension >1');
end;
if length(fn)>8;
    error('fn must be 8 or fewer chars');
end;
if ~isletter(fn(1));
    error('fn must begin with a letter');
end;

% suffx
if ~isstr(suffx);
    error('suffx not a string');
end;
if size(suffx,1)~=1;
    error('suffx has row dimension >1');
end;
if length(suffx)<2 | length(suffx)>3;
    error('suffx must be 2 or 3 chars');
end;
suffx=lower(suffx);
if ~any([strcmp(suffx,'rw')  strcmp(suffx,'eww')  strcmp(suffx,'lww')]);
    error(['suffx must be rw, eww or lww']);
end;


% path1
if ~isstr(path1);
    error('path1 not a string');
end;
if size(path1,1)~=1;
    error('path1 has row dimension >1');
end;


% Rev2011-04-04:  to handle path and file naming problems in linux
% % if ~isletter(path1(1));
% %     error('path1 must begin with a letter');
% % end;
% % if ~strcmp(path1(end),'\');
% %     error('path1 must end with a backslash');
% % end;
if ~(isletter(path1(1))  ||  strcmp(path1(1),'/'));
    error('path1 must begin with a letter or forward slash');
end;
if ~(strcmp(path1(end),'\')  ||  strcmp(path1(end),'/'));
    error('path1 must end with a backslash or forward slash');
end;



%-----  WRITE FILE

% Build path\filename
pf1 = [path1 fn '.' suffx];

% Open for writhe
fid1=fopen(pf1,'w');

% Write person, data and first year
fprintf(fid1,'%s\r\n',person);
fprintf(fid1,'%s\r\n',dwhen);
fprintf(fid1,'%s\r\n',int2str(yrfirst));

x=[round(x) ; -9999]; % append -9999 to measurements
g = num2str(x);
s = strjust(g,'left'); % left justify

[ms,ns]=size(s);
s=[repmat(blanks(1),ms,1) s  repmat(blanks(1),ms,1)];


for n = 1:ms;
    sthis=s(n,:);
    fprintf(fid1,'%s\r\n',sthis);
end;
fclose(fid1);
disp([pf1 ' sucessfully written']);
    





















    



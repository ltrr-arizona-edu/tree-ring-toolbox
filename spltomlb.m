function A=spltomlb()

% Convert a strung-out ascii vector back into its desired matrix form.

% D Meko 8-14-92

% Typical use:  Have a matrix in SPLUS;  Write it out of splus as
% a column vector in ascii;  Need to use the matrix in matlab.

% Missing values:  This routine assumes that all missing values are coded
% as a number equal to toobig, which is hard-coded below.  Note that SPLUS
% uses NA as missing value.  You need to edix the ascii file to convert
% all NA to toobig before calling this function.  This function then
% replaces all values of toobig-or-larger with NaN, the missing value
% code of matlab.  This song and dance needed because matlab cannot
% directly import the NA in a numeric matrix.

% The ascii file to be imported (name screen-prompted) must
% have a 3-character suffix -- as in DATA.DAT

toobig=1e32;  % anything as big or bigger is considered misssing data

d = input('NAME OF ASCII FILE TO INPUT: ','s');
eval(['load ',d]);
dl=length(d);     % Assume 4-char suffix (incl period), and cut it off.
d(dl-3:dl)=[];

n=input('DESIRED NUMBER OF COLS: ');

eval(['B=',d,';']);  % still a strung-out vector

%********     Put NaN for "missing" data

L=B>toobig-1;
b=NaN;
B(L)=b(ones(sum(L),1),:);


%*********   Reshape the strung-out vector into a matrix

eval(['m=length(',d,')/n']);
A=zeros(m,n);
A(:)=B;


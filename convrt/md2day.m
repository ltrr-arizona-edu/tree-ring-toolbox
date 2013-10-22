function day=md2day(m,d)
% md2day: matrix of month,day-of-month to sequential day numbers (366 in year)
% day=md2day(m,d);
% Last revised 12-10-00
%
% Written first for quick reference of which sequential day any given day-of-month is.
% 
%
%*** INPUT
%
% m (1 x ?)i month (1 to 12);
% d (1 x ?)i day of month
%
%
%*** OUTPUT
%
% day ( 1 x ?) sequential day (out of 366)
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%

% INITIALIZE REFERENCE MATRIX

n=[31 29 31 30 31 30 31 31 30 31 30 31]'; % n days in each month


% Check that m,d col vectors
[m1,n1]=size(m);
[m2,n2]=size(d);
if n1~=1 | n2~=1 | m1~=m2;
    error('m and d must be same size col vectors');
end;

% Check for invalid day and month numbers
if any(d<0) | any(d>31); 
    error('d must be between 1 and 31');
end;
if any(m<0) | any(m>12);
    error('m must be between 1 and 12');
end;


% Check for invalid day within 30-day month
for i =[4 6 9 11];
    L=m==i;
    if any(L);
        x=d(L);
        if any(x>30);
            error('31 days in a 30-day month');
        end;
    end;
end;
L=m==2; % feb
if any(L);
    x=d(L);
    if any(x>29);
        error('30 days in a February');
    end;
end;


% Ref for first days of months
c=cumsum(n);
n1=[1; c+1];
n1(13)=[];


day=(n1(m)+d-1);


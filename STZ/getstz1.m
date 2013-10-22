function [D,yrD]=getstz1(pf1,dtype,names)
% getstz1:  Get matrix of time series from stz-format storage
% [D,yrD]=getstz1(pf1,dtype,names);
% Last revised 11-26-01
%
% Get matrix of time series from stz-format storage.  Extracts time series matrix of any
% of several data types (e.g., ring width, core index
%
%*** INPUT
%
% pf1 (1 x ?)s name of input file
% dtype (1 x ?)s  code for data type
%   X = ring width
%   IX -- core index, standard
%   EX -- core index, residual
%   IT -- tree index, standard
%   ET -- tree index, residual
%   ZI -- site index, standard
%   ZE -- site index, residual
% names {} serie identifies (e.g., BOL18A) (see notes)
%
%*** OUTPUT
%
% D (mD x nD)r time serie matrix, NaN padded as needed
% yrD (mD x 1)i year vector
%
%*** REFRENCES -- NONE
%*** UW FUNCTIONS CALLE -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% names{}.  If data type ZI or ZE, names has length 1 and names{1}==[], as there is only
% one series of the desired type
%


% Load file
eval(['load ' pf1 ';' ]);
nmsc = cellstr(nms); % nms from char mtx to cell



switch dtype;
    
% store series names, data vector and years pointer    
case 'X'; % ring width
    nser=length(names);
    id=nmsc;
    x=X; 
    YRS = yrs;
otherwise;
    error('Yet to code options other than ring width in getstz1');
end;


switch dtype;
case 'ZI';
case 'ZE';
otherwise;
    yr1 =repmat(NaN,nser,1);
    yr2 = repmat(NaN,nser,1);
    for n = 1:nser;
        idthis = names{n};
        i1 = strmatch(idthis,id);
        if length(i1)~=1;
            error([idthis ' not there']);
        else;
            YRSthis = YRS(i1,:);
            yr1(n)=YRSthis(1);
            yr2(n)=YRSthis(2);
        end;
    end;
    yron = min(yr1);
    yroff=max(yr2);
    yrvec = (yron:yroff)';
    nyear = length(yrvec);
    yrD=yrvec;
    D=repmat(NaN,nyear,nser);
end;



% Pull and store series

switch dtype;
case 'ZI';
case 'ZE';
otherwise;
    for n = 1:nser;
        idthis = names{n};
        i1 = strmatch(idthis,id);
        YRSthis = YRS(i1,:);
        nyr = YRSthis(2)-YRSthis(1)+1;
        igo =YRSthis(3);
        isp =igo + nyr -1;
        xthis = x(igo:isp);
        jgo = YRSthis(1)-yron+1;
        jsp = YRSthis(2)-yron+1; 
        D(jgo:jsp,n)=xthis;
     
    end;
end;






function [x,yr]=getone1(pfnm,varnm,seqno);
% getone1: get one ring-width or index series from a .mat tree-ring storage file
% [x,yrx]=getone1(pfnm,varnm,seqno);
% Last revised 2-10-00
% 
% Gets ring-width or index data for a single core or tree from a .mat tree-ring
% storage file
%
%*** INPUT
%
% pfnm (1 x ?)s  path/filename of .mat storage file with input data
% varnm (1 x ?)s   type of variable:
%    X=ring width
%    IX = standard core index
%    EX = residual core index
%    IT = standard tree index
%    ET = residual tree index
% seqno (1 x 1)i sequence number series within storage file
%
%
%*** OUTPUT
%
% x (nx x 1)r  the desired time series
% yrx (nx x 1)i  year vector for x
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES -- none

varallow={'X','IX','EX','IT','ET'};
L=strcmp({varnm},varallow);

if ~any(L);
   error([pfnm ' does not contain ' varnm ]);
end;

if L(1) | L(2) | L(3); % if ring width or core index
   yrsnm = 'yrs'; % yrs holds the pointer data
else;
   yrsnm = [varnm 'yrs'];
end;

% Load the file
eval(['load ' pfnm ' ' varnm ' ' yrsnm]);

% Put the variables in specific file
eval(['x = ' varnm ' ;']);
eval(['yrx =  ' yrsnm ' ;']);


% Get the desired serie
n=seqno;
yr=(yrs(n,1):yrs(n,2))';
nyr=length(yr);
igo = yrs(n,3);
isp = igo+nyr-1;
x  = X(igo:isp);







%
%
%    
function Y=wy2mon01(W)
% wy2mon01:  monthly Calif DWR water year precip spreadsheet to .mat file
% wy2mon01;
% Last revised 6/25/01
%
% Convert a 13-col matrix of monthly data from water year to calendar year format (see notes)
% 
%*** INPUT
%
% W(? x 13)r year and data values for Oct_t-1 thru Sept_1
% 
%
%*** OUTPUT
%
% X (? x 13)r  year and data values for Jan - Dec, NaN filled
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% First use on Sacramento DWR to get N Sierras 8-station series in calendar year


[mW,nW]=size(W);
if nW~=13;
    error('W Must have 13 cols');
end;

% Water year data
yrW=W(:,1);



% Flesh out W to continuous years
yrgo = W(1,1);
yrsp = W(mW,1);
yrV = (yrgo:yrsp)';
nyrV = length(yrV);
V=repmat(NaN,nyrV,13);
V(:,1)=yrV;
irow = yrW-yrgo+1;
V(irow,:)=W;

% Trim all except monthly cols from V
V = V(:,2:13);

% Reorganize V into calendar year as X
B = V(:,1:3); % oct-dec
A = V(:,4:12); % jan-sept

yrX =[(yrV(1)-1); yrV];
A1=[(repmat(NaN,1,9)); A];
B1=[B; (repmat(NaN,1,3))];
X=[A1 B1];
[mX,nX]=size(X);
if all(isnan(X(1,:)));
    yrX(1)=[];
    X(1,:)=[];
    [mX,nX]=size(X);
end;
if all(isnan(X(mX,:)));
    yrX(mX)=[];
    X(mX,:)=[];
    [mX,nX]=size(X);
end;
Y=[yrX X];

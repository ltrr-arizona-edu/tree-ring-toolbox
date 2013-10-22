function m_elev(varargin);
% M_ELEV Contour elevation onto a map using a 1-degree database
%        M_ELEV contours elevations at 1000m intervals for the map.
%        M_ELEV(OPTN (,LEVELS) (,ARGS,...) ) lets you change various options.
%        if OPTN=='contour', contour lines are drawn. for OPTN=='contourf',
%        filled contours are drawn. LEVELS are the levels used, and ARGS
%        are optional arguments of line types, colors, etc. 
%
%        See also M_PROJ, M_GRID, M_COAST

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% 20/5/97 - moved registration over by 1/2 degree (seems to fit better)

global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;


if nargin==0,
 levels=[-7000:1000:-1000 1000:1000:5000];
 optn='contour';
 n_opt=1;
else
 if isstr(varargin{1}),
   optn=varargin{1};
 end;
 if nargin==1,
   levels=[-7000:1000:-1000 1000:1000:5000];
   n_opt=2;
 else
   if isstr(varargin{2}),
     levels=[-7000:1000:-1000 1000:1000:5000];
     n_opt=2;
  else
     levels=varargin{2};
     n_opt=3;
   end;
 end;
end;

load topo

blat=max(floor(MAP_VAR_LIST.lats(1)+.5),-89)-.5;
tlat=min(ceil(MAP_VAR_LIST.lats(2)+.5),90)-.5;
llong=floor(MAP_VAR_LIST.longs(1)+.5)-.5;
rlong=ceil(MAP_VAR_LIST.longs(2)+.5)-.5;

if rlong>360, rlong=rlong-360; llong=llong-360; end;
if llong<-360, rlong=rlong+360; llong=llong+360; end;

lts=(blat:tlat);
lgs=(llong:rlong);

if rlong<0,
  topo=topo(lts+90.5,lgs+360.5);
elseif llong<0 & rlong>=0,
  topo=topo(lts+90.5,[(360.5+llong:end) (1:rlong+0.5)]);
else
  topo=topo(lts+90.5,lgs+.5);
end;

[lg,lt]=meshgrid(lgs,lts);
[X,Y]=m_ll2xy(lg,lt,'clip','on');

hold on;
switch optn,
 case 'contour',
    contour(X,Y,topo,levels,varargin{n_opt:end});
 case 'contourf',
    i=isnan(X);
    topo(i)=NaN;
    [X,Y]=m_ll2xy(lg,lt,'clip','off');
    contourf(X,Y,topo,levels,varargin{n_opt:end});
end;  


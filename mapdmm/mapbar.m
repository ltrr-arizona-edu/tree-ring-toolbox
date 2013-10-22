function [lat1,lon1]=mapbar (lat0,lon0,rng,ndivs,units)
% mapbar: build map scale bar
% CALL: [lon1,lon1]=mapbar(lat0,lon0,rng,ndivs,units);
%
% Meko 9-7-97
%
%****************  IN 
%
% lat0(1 x 1)r  leftmost point (in longitude degrees) of bar
% 

rngorig = rng;

if strcmp(units,'km');
   rng=km2deg(rng);
elseif strcmp(units,'sm');
   rng=sm2deg(rng);
end

[lat1,lon1]=reckon(lat0,lon0,rng,90);
lat = [lat0;lat1];
lon= [lon0;lon1];

% Compute thickness of bar
x1 = getm(gca,'MapLatLimit');
x2= diff(x1)/60;

% Compute delta longitude
ydelta = diff(lon)/ndivs;

% Compute delta latitude
xdelta = diff(lat)/ndivs;

lonll=zeros(ndivs,1);
latll=zeros(ndivs,1);

% Compute lower left point for each division of bar
for n = 1:ndivs;
   lonll(n) = lon0 + (n-1)*ydelta;
   latll(n) = lat0 + (n-1)*xdelta;
   
end

lonul = lonll;
lonlr = lonll + ydelta;
lonur = lonlr;

latul = latll+x2;
latlr = latll + xdelta;
latur = latlr+x2;

% Plot bars
for n=1:ndivs;
   lat = [latll(n);latul(n);latur(n);latlr(n)];
   lon=[lonll(n);lonul(n);lonur(n);lonlr(n)];
   if rem(n,2)==0;
      hbar=fillm(lat,lon,[.8 .8 .8]);
   else
      hbar=fillm(lat,lon,[.3 .3 .3]);
   end
end


% Add text
%t1a=textm(latul(1)+x2/2,lonul(1),...
  % '0','Fontsize',8,'HorizontalAlignment','right',...
  % 'VerticalAlignment','baseline');
t1b=textm(latlr(ndivs),lonlr(ndivs)+x2/4,[num2str(rngorig) ' ' units],...
   'FontSize',8,...
   'HorizontalAlignment','left',...
   'VerticalAlignment','baseline');
   

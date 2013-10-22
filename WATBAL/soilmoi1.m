function meat=soilmoi1(p,pe,awcs,awcu,ssgo,sugo)
% soilmoi1:  end of month soil moisture from soil moisture accounting
% CALL: meat=soilmoi1(p,pe,awcs,awcu,ssgo,sugo);
%
% Meko  6-4-97
%
%*********************** IN *************************************
%
% p(nmos x 1)r monthly ppt (inches), assumed to begin with Jan of
%		first year and end with Dec of last year
% pe (nmos x 1)r monthly potential evapotranspiration (in); same
%		structure as p
% awcs (1 x 1)r available water capacity (in) in surface layer
% awcu (1 x 1)r available water capacity (in) in underlying layer
% ssgo (1 x 1)r starting soil moisture (in) in surface layer
% sugo (1 x 1)r starting soil moisture (in) in underlying layer
%
%
%*************************  OUT  *******************************
% meat{} -- all nmos x1 vectors
%  1- dels, soil moisture change in sfc layer
%  2- delu, soil moisture change in undrlying layer
%  3- del , soil moisture change in both layers
%
%  4- ss1, starting soil moisture, sfc
%  5- su1, starting soil moisture, underlying layer
%  6- s1, starting soil moisture, combined layers
%
%  7- ss2, ending soil moisture, sfc
%  8- su2, ending soil moisture, under
%  9- s2,  ending soil moisture, combined layers
% 10- smean, mean soil moisture, combined layers, for month  == (s1 + s2)/2
%
% 11-r, recharge, combined layers
% 12-pr, potential recharge
% 13-ro, runoff
% 14-pro, potential runoff
% 15-loss, loss
% 16-ploss, potential loss
%
% 17- et, estimated actual evapotranspiratin
%
%
%**************************  NOTES *****************************************
%
% Source: Palmer, Wayne C., 1965.  Meteorological Drought; Research Paper No. 45.
% US Dept of Commersce, Washington, DC.  Coded from equations and method described
% on pages 6-11.
%
%---------------
% One of several functions written with aim of getting easily modified Palmer Drought
% Index.  As is, the sequence of functions, from calling on down is:
%
% pdsicalc.m   -- hydacc1.m -- soilmoi.m
%
%----------------------------------
% pdsicalc.m computes the drought index.  
%
%----------------------------------------
% hydacc1.m controls some features of the
% hydrologic accounting used in the Palmer calculations.  For example, I have
% planned an option to allow modifying the monthly distrib of 'effective' ppt 
% depending on whether ppt is estimated to have fallen as rain or snow.  
% hydacc1.m also allows a pass-1 run to automatically come up with reasonable
% starting values for soil moisture for the first jan of the accounting period
%
% soilmoi1.m carries out the actual monthly accounting of soil moisture and other
% variables.



%********************* Check inputs


%----------------------  p, pe

[nmos,ntemp]=size(p);
if ntemp~=1; 
	error('p must be col vector');
end

[mtemp,ntemp]=size(pe);
if ntemp~=1; 
	error('pe must be col vector');
end
if mtemp~=nmos;
	error('pe must be same length as p');
end

if rem(nmos,12)~=0;
	error('length of p and pe  must be multiple of 12');
end

% ------------- awcs,awcu,ssgo,sugo

[mtemp,ntemp]=size(awcs);
if mtemp~=1 | ntemp~=1;
	error('awcs must be scalar');
end
[mtemp,ntemp]=size(awcu);
if mtemp~=1 | ntemp~=1;
	error('awcu must be scalar');
end
[mtemp,ntemp]=size(ssgo);
if mtemp~=1 | ntemp~=1;
	error('ssgo must be scalar');
end
[mtemp,ntemp]=size(sugo);
if mtemp~=1 | ntemp~=1;
	error('sugo must be scalar');
end

if ssgo>awcs | sugo>awcu;
	error('Starting soil water greater than awcs or awcu');
end

if awcu <=0 | awcs <=0;
	error('Avail water capacities must be greater than zero');
end

if ssgo<0 | sugo<0;
	error('Starting soil moisture levels must be non negative');
end


%******************** Size outputs

anan=NaN;
ss1=repmat(anan,nmos,1);
ss2=repmat(anan,nmos,1);
su1=repmat(anan,nmos,1);
su2=repmat(anan,nmos,1);

% Recharge
rs=repmat(anan,nmos,1);
ru=repmat(anan,nmos,1);

% Runoff
ro=repmat(anan,nmos,1);

% Net loss to Evapotranspiration
es=repmat(anan,nmos,1);
eu=repmat(anan,nmos,1);

% Change in soil moisture
dels=repmat(anan,nmos,1);
delu=repmat(anan,nmos,1);

%**************************** Water balance

d=pe-p; % deficit=excess of potential evapo over precip
awc=awcu+awcs; % combined water capacity of two layers

%************************  Loop over Months

ss1this=ssgo; % initialize starting soil moisture, sfc, first month
su1this=sugo; % .... underlying

for n = 1:nmos;
	
	dthis = d(n); % pe-p for this month
	
	% Sfc Layer
	sempty = awcs-ss1this;  % how much could sfc layer take in
	if dthis>=0; % if pe exceeds precip
		dels(n)=-dthis; % tentatively set change in soil moisture to pe-p
		if dthis>ss1this; % if pe-p for the month exeeds avail water in sfc layer
			dels(n) = -ss1this; % sfc layer can only lose what it has
		end
		rs(n)=0; % no recharge to sfc layer
		ro(n)=0; % no runoff this month
		% Net  Loss from sfc layer
		if dels(n)<0;
			es(n) = -dels(n);
		else;
			es(n) = 0;
      end
      excess=0;
	else; % ppt exceeds pe
		dels(n)= min(sempty,-dthis);
		rs(n)=dels(n);
		excess=(-dthis) - dels(n);
		es(n)=0;
	end
	ss1(n)=ss1this;
	ss2(n)=ss1this+dels(n);
	ss1this = ss2(n);


	% Underlying layer

	uempty = awcu-su1this;  % how much could sfc layer take in
	if excess<=0;  % No excess from upper layer
		eu(n) = (dthis - es(n)) * (su1this/awc); % 'loss' from underlying layer
      eu(n)=min(eu(n),su1this);
      if eu(n)<0;
         eu(n)=0;
      end
      ru(n)=0;
		ro(n)=0;
		delu(n) = -eu(n); % change in underlying soil moisture
		
	else; % some excess from overlying layer
		eu(n)=0; % no 'loss' from underlying layer
		delu(n)=min(uempty,excess);
		ru(n) = delu(n);
		if excess>uempty;
			ro(n)=excess-uempty;
		else;
			ro(n)=0;
		end

	end
	su1(n)=su1this;
	su2(n)=su1this+delu(n);
	su1this = su2(n);
end

%***********  Vector math for other quantities

del=delu+dels; % change in soil moisture, combined layers
et = p - ro - del;  % evapotranspiration
r = rs + ru; % recharge
loss = es + eu;
s1 = ss1+su1;  % start of month soil moisture in both layers
s2=  ss2+su2; % end of month soil moisture in both layers
smean=(s1+s2)/2; % mean of starting and ending soil moisture


%******************** Potential recharge, loss and runoff
pr = awc - s1; % potential recharge

% Potential loss
dope1=[pe ss1]';
plosss=(min(dope1))'; % pot loss, sfc layer
plossu= (pe-plosss) .* (su1/awc);
ploss=plosss+plossu;

pro = awc - pr; % potential runoff

%******************  STORE RESULTS IN CELL VARIABLE TO PASS BACK TO CALLING FUNCTION

meat{1}=dels;
meat{2}=delu;
meat{3}=del;

meat{4}=ss1;
meat{5}=su1;
meat{6}=s1;

meat{7}=ss2;
meat{8}=su2;
meat{9}=s2;
meat{10}=smean;

meat{11}=r;
meat{12}=pr;

meat{13}=ro;
meat{14}=pro;

meat{15}=loss;
meat{16}=ploss;

meat{17}=et;


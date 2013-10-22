function e=pltrec3(yr4,modnum,EV,RE,NV,ktrans)
% e=pltrec3(yr4,modnum,EV,RE,NV,ktrans)
% See spatrec1.m for instructions

a=NaN;

% Compute time series of EV,RE,n1,n2,n3
for n = 1:length(modnum);
	L1=modnum(n);
	zev(n) = EV(modnum(n));
	zre(n) = RE(modnum(n));
	zn1(n) = NV(modnum(n),3); % available trees
	zn2(n) = NV(modnum(n),4); % potential pc-amp predictors
	zn3(n) = NV(modnum(n),5); % final pc-amp predictors
end

% Check current figure background color and change if necessary
curcolor=get(gcf,'Color');
black=all(curcolor==0);
if ktrans==0 & black==1; % OK for black bkgd slide
elseif ktrans==0 & black==0; % change background
	whitebg;
elseif ktrans==1 & black==0; % OK for white background transparency
elseif ktrans==1 & black==1; % change background
	whitebg
else
	error('unknown combination of ktrans and black');
end


subplot(2,1,1)
h1=plot(yr4,zn1,':',yr4,zn2,'--',yr4,zn3,'-');
set(h1,'Linewidth',1)
%set(hh1,'Linestyle',':','Linewidth',1)
ylabel('Number of Variables')

% Change line colors, thicknesses depending on ktrans
if ktrans==0; % want black bkg figure.  Colors OK as is
elseif ktrans==1; % want white bkg figure
	set(h1(3),'Linewidth',2,'Color','k')
	set(h1(2),'LineStyle','-','Linewidth',2,'Color','m')
	set(h1(1),'LineStyle','-','Linewidth',2,'Color','r')
else
	error('ktrans must be 0 or 1');
end
legend('Trees','Potential Predictors','Selected Predictors')
title('VARIATION OF PREDICTOR MAKEUP WITH TIME')

%********************************


subplot(2,1,2)
h2=plot(yr4,zev,'-',yr4,zre,'--')
set(h2,'LineWidth',1)
ylabel('Accuracy')
xlabel('Year')

% Change line colors, thicknesses depending on ktrans
if ktrans==0; % want black bkg figure.  Colors OK as is
elseif ktrans==1; % want white bkg figure
	set(h2(1),'LineStyle','-','Linewidth',2,'Color','k')
	set(h2(2),'LineStyle','-','Linewidth',2,'Color','m')
else
	error('ktrans must be 0 or 1');
end


title('VARIATION OF RECONSTRUCTION ACCURACY WITH TIME')
legend('R-Squared','RE')




function recplot(yhat,se,ymean,yr,k5)

% Plot segment of time series, optionally with error bars

%*****   INPUT ARGS
%
%  yhat (m1 x 1)  the time series, a regression predicted series
%  se  (1 x 1)  std error of the estimate
%  ymean (1 x 1)  calib-period mean
%  yr(1 x 2) first, last years of vector yhat
%  k5 (1 x 1)  "Y" or "N"   -- keyboard mode pause on or off

ymall=mean(yhat);
years=(yr(1):yr(2))';

k3=[];
k3=input(['FIRST PLOT YEAR ',int2str(yr(1)),'? Y/N  [Y] :'],'s');
if isempty(k3), k3='Y';  end;
if k3=='Y'
	 yrplot(1)=yr(1);
else
	yrplot(1)=input('FIRST YEAR FOR PLOT? ');
end

k3=[];
k3=input(['LAST PLOT YEAR ',int2str(yr(2)),'? Y/N  [Y] :'],'s');
if isempty(k3), k3='Y';  end;
if k3=='Y'
	 yrplot(2)=yr(2);
else
	yrplot(2)=input('LAST YEAR FOR PLOT? ');
end

L1 = [yrplot(2)<=yrplot(1)  yrplot(1)<yr(1)  yrplot(2)>yr(2)];
if any(L1)
	clc; home;
	disp('Aborting recplot.m: start,end year for plot impossible!');
	pause
	return
end


yrp = (yrplot(1):yrplot(2))';
L1=  years >= yrplot(1)   & years <= yrplot(2);
y=yhat(L1);


k2=1;  % control loop for varying std err bars
while k2 ~=7;  %stay in loop
k2=menu('ERROR-BAR OPTION','+- TWO SE',...
'+- ONE SE','NO SE BARS','USER-DEFINED',...
'NO MORE PLOTTING');

if k2==1;  % plot bars at 2 se around estimate
	fact=se*2.0;
elseif k2==2;  % one se bars
	fact=se;
elseif k2==3; % no bars
	fact=0;
elseif k2==4; % user-set multiple of se 
	ff=input('MULTIPLE OF SE FOR ERROR BARS:  ');
	fact=ff * se;
elseif k2==5;
	 k2=7;
	break
end;




k1=1;
while k1~=5;  % stay in loop
clg; clc;
k1=menu('MEAN-LINE OPTION','FULL-LENGTH','CALIB PERIOD',...
'PLOT SUBPERIOD','OTHER','QUIT');

if k1==1,
	ym= ymall; 
elseif k1==2,
	ym=ymean(1,1); 
elseif k1==3,
	ym = mean(y); 
elseif k1==4,
	mnper=input('PERIOD FOR MEAN: ');
	L2 = years >= mnper(1) & years <= mnper(2);
	ym = mean(yhat(L2));
elseif k1==5,
	
	break
end

n=length(y);
ybar=ym(ones(n,1),:);  % dupe mean into a cv


if k2==3, 
	symb='--i';
else
	symb='--g';
end

biggy=max([max(y)  max(ybar+fact)]);
small=min([min(y)  min(ybar-fact)]);

rng=abs(biggy-small);
pymin=small - (rng/20);
pymax=biggy + (rng/20);

plot(yrp,y,'-',yrp,ybar+fact,symb,yrp,ybar-fact,symb,...
yrp,ybar,'-g',yrp(1),pymin,'i',yrp(1),pymax,'i');
xlabel('YEAR')

if k5=='Y', keyboard;
else, pause, end;
clg


end;  % of while on k1
end; % of while on k2

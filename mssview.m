function nmiss=mssview
%
% Graphic picture of where monthly pcp data is missing for a station
%
%
[file1,path1]=uigetfile('*.mat','File with monthly climatic data'); 
ff=[path1 file1];
eval(['load ' ff]);


% assume yr and data in X, Y, Z or Y1
if ~exist('Y1') & ~exist('Z') & ~exist('X')& ~exist('Y');
	error('Y1,X,Z or Y not in this mat file')
end

%-- Check that only one of the allowable climate variable names is used
Lhere = [exist('Z') exist('Y1') exist('X') exist('Y')];
if sum(Lhere)~=1;
	error('More than one of Z,Y1,X or Y is in the loaded file');
end


%----  Put the monthly climate data in a matrix called Z
if exist('Y1');
	Z=Y1;
end
if exist('X');
	Z=X;
end
if exist('Y');
   Z=Y;
end

if exist('Z');
	% no action needed, Z is what you want
end

[m1,n1]=size(Z); % Z has m1 years, n1 columns

% Put data part in A, year in yr
A=Z(:,2:13);
yr = Z(:,1);

% Any NaN's in A?
L1=isnan(A);

if any(any(L1)); % some missing values in A
   nmiss=sum(sum(L1));
else;
   nmiss=0;
end;


figure(1)
clf
hold on

txt1=['Missing Data for Station: ' strtok(file1,'.')];
txt2=[int2str(yr(1)) '-' int2str(yr(m1))];


% Make missing values ==1 for Jan, 2 for feb, etc
for n=1:12; % loop over months;
	L2=L1(:,n);
	if any(L2);
		i=find(L2);
		t=yr(i);
		y=L2(i)*n;
      plot(t,y,'rx');
      
	else
	end
end

grid
title([txt1 ' ' txt2]);

ylabel('Month')
xlabel('Year')
figure(1)
set(gca,'Xlim',[yr(1) yr(m1)],'Ylim',[0 13]);
if length(yr)<10;
   set(gca,'XTick',[min(yr):max(yr)]);
end

zoom xon
hold off

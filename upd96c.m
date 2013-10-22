function upd96b
%
% Update monthly pcp or tmp file thru Oct 1996 by
% promping for Sept-Oct 96 data


clear

[file,path]=uigetfile('*.mat','File to be updated')
pf=[path file];
eval(['load ' pf]);


% Usually, my monthly data is in matrix Z, but sometimes in Y1
if exist('Y1') & ~exist('Z')
	Z=Y1;
end


[m1,n1]=size(Z);
clc
years=Z(:,1);
yr1=Z(1,1);
yr2=Z(m1,1);
disp(['First, last years ',int2str(yr1),' ',int2str(yr2)]);
disp(' ');

a=NaN;
x=a(:,ones(13,1));


% Find current 1996 data, thru Aug, and put in slots in x
i1=find(years==1996);
if isempty(i1);
	error('No current 1996 data in file');
end

% Put current 1996 data in x
x=Z(i1,:);


k=menu('Choose One',...
'Manually key in Sept-Oct 96 data',...
'Add NaN years to bring matrix through 1996');
if k==1;
	x9=input('Sept value: ');
	x10=input('Oct value: ');

	disp('')
	x(10:11)=[x9 x10];
	disp('Here is the 1996 data you keyed in')
	disp('')
	fprintf('%4.0f',x(1));
	for n=10:11;
		fprintf('%6.2f',x(n));
	end
	fprintf('\n');
elseif k==2;
	error('just kidding -- do not choose this option in upd96.m');
end
		
Z(i1,:)=x;
eval(['save ' pf  ' Z']);




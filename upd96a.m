function update96
%
% Update monthly pcp or tmp file thru Aug 1996 by
% promping for Jan-Aug 96 data


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
yr1=Z(1,1);
yr2=Z(m1,1);
disp(['First, last years ',int2str(yr1),' ',int2str(yr2)]);
disp(' ');

a=NaN;
x=a(:,ones(13,1));


k=menu('Choose One',...
'Manually key in Jan-Aug 96 data',...
'Add NaN years to bring matrix through 1996')
if k==1;
	x(1)=1996;
	x1=input('January value: ');
	x2=input('February value: ');
	x3=input('March value: ');
	x4=input('April value: ');
	x5=input('May value: ');
	x6=input('June value: ');
	x7=input('July value: ');
	x8=input('August value: ');



	disp('')
	x(2:9)=[x1 x2 x3 x4 x5 x6 x7 x8];
	disp('Here is the 1996 data you keyed in')
	disp('')
	fprintf('%4.0f',x(1));
	for n=2:9;
		fprintf('%6.2f',x(n));
	end
	fprintf('\n');
elseif k==2;
	need = 1996-yr2;
	if need==0;
		error('Series already has a 1996 line')
	end
	x = x(ones(need,1),:);
	yrvect=    ((yr2+1):1996)';
	x(:,1)=yrvect;

end
		
Z=[Z;x];
eval(['save ' pf  ' Z']);




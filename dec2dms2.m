function Y=dec2dms2(X,kopt)
% long-lat in decimal degrees to lat-long file in 
% (deg, min, sec) 
% Works correctly only for regions in N Hemispher, West longitude
% By D. Meko, 4-30-96

%*******   INPUT ARGS ****************
%
% X (mX,2) - dec degrees of longitude (col 1) and
%		latitude (col 2) for points
% kopt (1 x 2)i
% 1  ==1: west longitude negative in input
%	  ==2: west long. not set negative in input
% 2  ==1: some user dialog
%	  ==2: quickie, no user dialog
%
%************* OUT ARGS *************
%
% Y (mY x 6)r - (deg, min,sec) of latitude, followed by longitude 
%
% Optional output -- ascii file with converted long,lats.
% User is prompted for file name
%

%warndlg('Works only in NH, W long')


[mX,nX]=size(X);
if nX ~=2
	error ('X must be 2-column matrix')
end


[m1,n1]=size(kopt);
if m1~=1 | n1~=2
	error('kopt must be 1 x 2')
end

if kopt(1)==1; % input decimal longitude negative for west; change sign
	L=X(:,1)>0;
	if any(L),
		error('All longitude input were not negative')
	end
	X(:,1)=-1.0 * X(:,1);
elseif kopt(1)==2,
	% no action
else
	error('kopt(1) must be 1 or 2')
end



a=NaN;
Y = a(ones(mX,1),ones(6,1));

Y(:,1)=floor(X(:,2)); % degrees latitude
d1=(X(:,2)-Y(:,1))*60;
Y(:,2) = floor(d1); % minutes latitude
d2 = (d1 - Y(:,2))* 60;
Y(:,3)=round(d2);  % seconds latitude

Y(:,4)=floor(X(:,1)); % degrees longitude
d1=(X(:,1)-Y(:,4))*60;
Y(:,5) = floor(d1); % minutes longitude
d2 = (d1 - Y(:,5))* 60;
Y(:,6)=round(d2);  % seconds longitude


if any(any(Y)<0),
	error('An output Y value is negative')
end



if kopt(2)==1; % dialog- mode
	k1=input('Want an ascii output file? [Y]/N', 's')
	if isempty(k1) | k1=='y' | k1=='Y'
		file1=uiputfile('*.dat','Desired output ascii long-lat file: ')
		fid1=fopen(file1,'wt')
		for n=1:mX;
			z=[Y(n,1) Y(n,2) Y(n,3) Y(n,4) Y(n,5) Y(n,6)]';
			fprintf(fid1,'%5d %5d %5d %5d %5d %5d\n',z);
		end
		fclose(fid1)
	elseif k1=='N' | k1=='n'
		% no action
	else
		error('Response should be Y or N')
	end
end

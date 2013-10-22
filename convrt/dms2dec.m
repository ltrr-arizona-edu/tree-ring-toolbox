function Y=dms2dec(X)
% dms2dec: lat/long in 'dms' to long/lat in decimal degrees
% CALL: Y=dms2dec(X);
%
% Meko, 5-8-94
%
%*******   INPUT ARGS ****************
%
% X (mX,6 or 4) - deg, min,(sec) of latitude, followed by longitude 
%		Longitude degrees might be coded as negative if W longitude,
%		or might not be.  Program will prompt for information on this
%		coding.  
%
%************* OUT ARGS *************
%
% Y (mX,2) - decimal degrees of longitude, latitude for points specified
%   by X.  Order of longitude, latitude has been reversed from X so that
%   "x-coordinate" column is now before "y-coordinate" column.  Decimal 
%   longitude in Y is expressed as negative number, because autocad and
%   other mapping coord systems have x "increasing" to right (East).
%
% Optional output -- ascii file with converted long,lats.
% User is prompted for file name
%
%
%*****************  NOTES *******************************
%
% Works correctly only for regions in N Hemispher, West longitude
%
%******************  END OF OPENING COMMENTS ********************



[mX,nX]=size(X);

warndlg('Works only in NH, W long')

% Check whether input file contains seconds; store varibles
if nX==6; % degrees, minutes, seconds in input file
	LTd =X(:,1); % latitude degrees
	LTm =X(:,2);
	LTs =X(:,3);
	LGd =X(:,4); % longitude degrees
	LGm =X(:,5);
	LGs =X(:,6);
	X1 = X(:,[1 2 3 5 6]);
	% degrees lat, and minutes, seconds should never be coded negative
	if any(any(X1<0)),
		error('Negative value found in minutes or seconds')
	end
	% Check that minutes and seconds within range
	if any(any(X(:,[2 3 5 6])>60)),
		error('A minutes or seconds value greater than 60')
	end
elseif nX==4 % degrees, minutes only in input file
	LTd =X(:,1); % latitude degrees
	LTm =X(:,2);	
	LGd =X(:,3); % longitude degrees
	LGm =X(:,4);
	X1=X(:,[1 2 4])
	% lat/deg or lat,long mins, seconds should never be coded negative
	if any(any(X1<0)),
		error('Negative value found in minutes')
	end
	% Check that minutes  within range
	if any(any(X1(:,[2 3])>60)),
		error('A minutes value greater than 60')
	end
else
	error('X must have either 4 cols or 6 cols')
end

% check that longitude not mistakenly put in first col of X
if any(abs(X(:,1))>90),
	error('An input latitude degrees is greater than 90')
end


% Change sign of Longitude degrees before conversion if the
% the sign is already negative for West longitude
L1 = all(LGd>0);
L2= all(LGd<0);
s = input('Is W long coded as negative in X? [Y]/N ','s')
if isempty(s) | s=='y' | s=='Y'
	s='Y'
elseif s=='n' | s=='N'
	s='N'
else
	error('Answer should be Y or N or y or n')
end

if s=='Y'
	if L2~=1,
		error('All input LongDeg values are not negative')
	end
	LGd=-1.0*LGd;
else
	if L1~=1,
		error('All input LongDeg values are not positive')
	end
end


if nX==6
	Y(:,1)= -(LGd+LGm/60+LGs/3600);
	Y(:,2)=  LTd+ LTm/60+LTs/3600;
else
	Y(:,1)= -(LGd+LGm/60);
	Y(:,2)=  LTd+LTm/60;
end



k1=input('Want an ascii output file of the x, y coords? [Y]/N', 's')
if isempty(k1) | k1=='y' | k1=='Y'
   [file1,path1]=uiputfile('*.dat','Desired output ascii long-lat file: ');
   pf1=[path1 file1];
	fid1=fopen(pf1,'wt')
	for n=1:mX;
		z=[Y(n,1) Y(n,2)]';
		fprintf(fid1,'%7.2f  %7.2f\n',z);
	end
	fclose(fid1)
elseif k1=='N' | k1=='n'
	% no action
else
	error('Response should be Y or N')
end

button1=questdlg('Want an ascii file coords and point # as 3rd column?');
switch button1
case 'No';
   % No action needed
case 'Cancel';
   % No action needed
case 'Yes';
   [file1,path1]=uiputfile('*.dat','Desired output ascii file name for x,y,n data: ');
   pf1=[path1 file1];
   fid1=fopen(pf1,'w');
   n = (1:mX)';
   z=[Y(:,1) Y(:,2) n];
   fmt1='%7.3f %7.3f %4.0f\n';
   fprintf(fid1,fmt1,z');
   fclose(fid1);
end   

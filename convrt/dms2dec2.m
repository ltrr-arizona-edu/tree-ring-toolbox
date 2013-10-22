function Y=dms2dec2(X)
% Y=dms2dec2(X)
% lat-long file in (deg, min, sec) or (deg,min) to decimal deg
% Works correctly only for regions in N Hemispher, West longitude
% Differs from dms2dec.m in that no user prompted info.  Written
% for automated run in making flow database for dettinger HCDN data
%
% By D. Meko, 10-17-96
%
%*******   INPUT ARGS ****************
%
% X (mX,6 or 4) - deg, min,(sec) of latitude, followed by longitude.
%		Longitude degrees might be coded as negative
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
[mX,nX]=size(X);


% Check whether input file contains seconds; store varibles
if nX==6; % degrees, minutes, seconds in input file
	LTd =X(:,1); % latitude degrees
	LTm =X(:,2);
	LTs =X(:,3);
	LGd =X(:,4); % longitude degrees
	LGm =X(:,5);
	LGs =X(:,6);
	X1 = X(:,[1 2 3 5 6]);
	% minutes, seconds should never be coded negative
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
if L2
	LGd=-1.0*LGd;
end


if nX==6
	Y(:,1)= -(LGd+LGm/60+LGs/3600);
	Y(:,2)=  LTd+ LTm/60+LTs/3600;
else
	Y(:,1)= -(LGd+LGm/60);
	Y(:,2)=  LTd+LTm/60;
end




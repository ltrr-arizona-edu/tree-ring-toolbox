function fpfmon2(Y,fn,tit,dtype)
% Print summary mean-ratio statistics
%
% D Meko 10-27-94
%
%********************  IN ARGS  ***************************
%
% Y (mY x 13)r   col 1 holds the "sequence number" of a predictor
%	station.  Cols 2-13 hold data values for jan-dec. Some rows
%  might be all zeros;  The function drops these off
% fn (?)s  name for output ascii file
% tit (1 x 64)s  title 
% dtype (1 x 1)i data type: 1=mean ratios, 2=correlations;
%		3=sample size (number years) for mean ratio;  Need this to
%		decide on %d vs %f formatting; and on truncating real data
%
%***************************  OUT ARGS -- NONE
%
%********************  NOTES  **********************
% 
% Mean ratios printed to 3 digits to right of decimal; Correlation
% to 2 digits.  Sample size printed as integer


% Check input mtx
[mY,nY]=size(Y);
if(nY~=13);
	error('ppt matrix should be 13 cols')
end


% Compute the number of predictor stations; lop off null rows;
L1=Y(:,1)>0;
ns = sum(L1);
Y=Y(L1,:);
[mY,nY]=size(Y);


% Pull out station sequence 
ID= round(Y(:,1));
Z = Y(:,2:13);
[mZ,nZ]=size(Z);


% If data are mean ratios  round to 3 digits right of decimal
% If data are correls, round to 2 decimal points to right
if dtype==1 
	Z = round(Z*1000)/1000;
elseif dtype==2
	Z = round(Z*100)/100;
end
	

% open file output file for writing
fid=fopen(fn,'w');

% Build title
blnk=blanks(64);  % make 1 x 64 string of blanks
lt=length(tit); % Number of characters in string argument for title
if lt > 64
	error('Make title 64 or fewer characters')
elseif lt==64
	tit1=tit;  % use title as is
else;  % center title
	ndiff=64-lt;
	nn=fix(ndiff/2);
	remmy = rem(ndiff,2);
	if remmy==0
		tit1 = [blanks(nn)  tit blanks(nn)];
	else
		tit1=  [blanks(nn+remmy) tit blanks(nn)];
	end
end
end


% Print title followed by blank line 
fprintf(fid,'%s',tit1);
fprintf(fid,'\n\n');


% Print header line after tailoring to dtype 
hdr1 = 'Stn    J      F      M      A      M      J';
hdr2 =  '      J      A      S      O      N      D';
if dtype==1 | dtype==2
	 hdr = [hdr1 hdr2];
else
	hdr='Stn   J   F   M   A   M   J   J   A   S   O   N   D';
end
fprintf(fid,'%1s',hdr);   
fprintf(fid,'\n\n');  % position at next line for data


for i=1:mY;    % loop over predictors
	fprintf(fid,'%3d',ID(i));  % print the ID number
	% Print the remaining values in a row
	z=Z(i,:);  % put this rows data in a vector
	if dtype==1 
		fprintf(fid,'%7.3f',z);
	elseif dtype ==2
		fprintf(fid,'%7.2f',z);
	elseif dtype==3
		fprintf(fid,'%4d',z);
	else 
		error('Hey, dtype is not equal to 1,2 or 3')
	end
	fprintf(fid,'\n'); % new line
end



fclose(fid);


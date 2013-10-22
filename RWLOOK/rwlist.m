function rwlist(X,Y)
%
% List common period of two time series in ascii file
% aatemp.dat in the current directory. Used by rwlook to allow
% easy viewing of the questionable series (X) and the known series
% (Y). Three-column list: year, X, Y
%
%******************  IN ARGS ********************************
%
% X (mX x 2) the first time series, year in col 1
% Y (mY x 2) the second ...

[mX,nX]=size(X);
[mY,nY]=size(Y);

tx=X(:,1);
ty=Y(:,1);
t1=max([tx(1) ty(1)]);
t2=min([tx(mX) ty(mY)]);

LX = tx >= t1 & tx<= t2;
LY = ty >= t1 & ty<= t2;

yr=(t1:t2)';
x=X(LX,2);
y=Y(LY,2);
j=(1:length(yr))';
xmean=mean(x);
ymean=mean(y);
y = y * xmean/ymean;  % scale master by ratio of means of sample to
%		master to ease comparison of listing of the two series

% Make shifted sequential numbering vector so that aatemp.dat
% shows the ring-count number of each ring in the undated
% series.  
if tx(1)>=ty(1); % if sample begins in same year or later
	% than master
	jj=j;  % ring "1" in parens of aatemp.txt is first year on
		% sample
else
	jj=j+(ty(1)-tx(1));
end


Z=[j jj yr x y];

fid=fopen('aatemp.dat','wt'); % open file for writing as text
fprintf(fid,'%6.0f %6.0f %5.0f %6.0f  %6.0f\n',Z');
fclose(fid);



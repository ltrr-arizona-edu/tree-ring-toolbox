function os3asci
%
% Put data from .mat output from screen3.m into and ascii
% file suitable for surfer to read to make maps of variables
% D Meko 12-16-96
%
% 

% should now have the following matrices in the function space.  Note
% that all cols are NaN except col 1
%
% col 1 --sequence
% col 2 -- x plotting coord
% col 3 -- y plotting coord
%
% NY (nT x 1)i number of years used to fit final OE model
% OB (nT x 1)i  order of B-operator in model
% OF (nT x 1)i  order of F-operator in model
% K2 (nT x 1)i  significant signal?  meaning are any of the
%		parameters of the full-period model signifantly different
% 		from zero (2+ sdevs); 1=yes 0=no
% G (nT x mxs)i number of significant (99% level) autocorrelation
%		coefficients  of residuals of model
% H (nT x 1)i number of significant (99% level) cross-
%		correlations between input u(t) and residuals e(t)
% R1 (nT x 1)r  corr coefficient between y(t) and u(t)
% R2 (nT x 1)r  corr coefficient between y(t) and [B/F] u(t)
% S (nT x 1)r variance-explained fraction for best OE model
%		Computed as 1 - (var(residuals)/var(y)) 
%		for final model fit to entire data period

% Get .mat file
[file1,path1]=uigetfile('*.mat','Output file from screen3.m');
pf1=[path1 file1];
eval(['load ' pf1]);

% Get row size
mr=size(H,1);

% Get x,y plot coordinates for tree sites
[file3,path3]=uigetfile('*.dat','x,y plot coordinates');
pf3=[path3 file3];
eval(['load ' pf3]);
eval(['C= ' strtok(file3,'.')]);
mC=size(C,1);
if mC ~=mr;
	error('x,y coord file must be same row-size as data file');
end

% File for ascii output
[file2,path2]=uiputfile('*.dat','Ascii out file');
pf2=[path2 file2];
fid2=fopen(pf2,'w');

% Combine first cols of matrices
X1=[(1:mr)' C  NY(:,1) OB(:,1) OF(:,1) K2(:,1) G(:,1) H(:,1)];
X2=[R1(:,1) R2(:,1) S(:,1)];
X=[X1 X2];

% format for a line of output ascii file
fmt2a='%4.0f %8.2f %6.2f %3.0f   %2.0f %2.0f  %1.0f    %2.0f %2.0f ';
fmt2b='%6.3f %6.3f %6.3f\n';
fmt2=[fmt2a fmt2b];

fprintf(fid2,fmt2,X');

fclose all

% spatout1.m -- script file for summary table (tab-delimited file)
% of single-site reconstructons
%
% D Meko 12-24-95
%
%*** U-W FUNCTIONS CALLED
%
% tabout2.m
% INSTRUCTIONS. 
% - run spatrec1.m, including crossvalidation and saving results
%	 to some .mat file
% - clear the workspace
% - load the .mat file
% - run spatout1.m
%

% Columns in output table
% 1 - site (or tree) number; sequential numbering that is same as in
%	 tree-index matrix columns
% 2-3 -- start and end valid years of single-site recons
% 4-7 -- estimated regression parameters, in following order
%		constant
%		t-1
%		t
%		t+1
% 8 -- regression R-squared
% 9 -- overall-F for regression equation
% 10 -- p-value for overall F
% 11 -- crossvalidation squared correl coef between observed and
%		predicted
% 12 -- crossvalidation reduction-of-error statistic

[mt,nt]=size(STATS);
s1a = Cwgt';
s1a = s1a(:,[1 3 2 4]);
s1=[(1:mt)' Pyrs s1a];
s2=STATS(:,[1 3 4]);
s3=[R2ss REss];

S=[s1 s2 s3];
hd1=['no ';'yr1';'yr2';'con';'t-1';' t ';'t+1'];
hd2=['R2 ';'F  ';'p-v';'R2v';'RE '];
tS=[hd1; hd2];
fmtd1='%2.0f\t%5.0f\t%5.0f\t%5.2f\t%5.2f\t%5.2f\t';
fmtd2='%5.2f\t%5.2f\t%5.2f\t%5.4f\t%5.2f\t%5.2f\n';
fmtd=[fmtd1 fmtd2];
fmth='%3s';
tabout2(S,tS,fmtd,fmth)


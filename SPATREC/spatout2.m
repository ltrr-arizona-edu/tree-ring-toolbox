% spatout2.m -- script file for summary table (tab-delimited file)
% of spatial reconstructions
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
% - run spatout2.m
%

% Columns in output table
% 1 - model number; squentially modeled, earliest first
% 2-3 -- start and end valid years of model
% 4-6-- predictor data-set reduction
%		4 - number of single-site predictors available
%		5 - number of PC amplitudes as potential predictors
%		6 - final number of PC amplitudes as predictors
% 7- regression R-squared
% 8 - validation squared r between actual and reconstructed
% 9 - validation reduction-of-error statistic
% 10 - rmse of prediction for validation data 

a = NaN;
S = a(ones(nmods,1),ones(10,1)); % allocate
s1 = [(1:nmods)'  P2 NV(:,3:5)];
s2 = [RSQ2 RSQ3  RE PEbar];
S=[s1 s2];

hd1=['no ';'yr1';'yr2';'NTR';'NPP';'NUP'];
hd2=['R2 ';'r2 ';'RE ';'RMS'];
tS=[hd1; hd2];
fmtd1='%2.0f\t%5.0f\t%5.0f\t%5.0f\t%5.0f\t%5.0f\t';
fmtd2='%5.2f\t%5.2f\t%5.2f\t%5.3f\n';
fmtd=[fmtd1 fmtd2];
fmth='%3s';
tabout2(S,tS,fmtd,fmth)


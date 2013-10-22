function screen7
% screen7:  tsm of chronologies selected in screen6.m; optional negation; ascii sov
%  Optional whitening by AR models; ascii long/lat plot file
% CALL:  screen7
%
% Meko 8-19-97
%
%*****************  IN *************************
% 
% tsm?.mat mother time series matrix available
% ?p#os6.mat screen6.m output file available
% crnxy?.dat ascii x,y long-lat file of mother set
% 
% Screen prompts:
%   start, end year for output time series matrix
%   whether to change sign of flagged series with negative clim-tree r
%   whether want AR residuals, and if so, max AR order to consider
%
%
%**************** OUT *************************
%
% .mat file with 
%   output time series matrix X, year in col 1
%   id string matrix N -- 
%   Ineg variables columns in output X of any series negated
%   output AR(1) or AR(2) residuals matrix same size as X; name E
%   AR modeling info in K1, W (see comments in text)
%   xycoords  plotting long, lat
%
% Ascii column vector of time series X, repeated col by col. One line
%  header with start year, end year , and 'original' or 'ar residuals'
%
%***************** NOTES ************************
%
% Sign changes, if any, are done before any AR modeling.  
% 
% AR modeling is done on full-set chronologies, after sign change (if any),
% not on the final output matrix period of data.  Of course, the full-set
% and output tsm might cover the same period.  Reasoning is to use as much
% data as possible to model dependence on past values.
%
% Residuals E will have 1 or more startup values equal to chron mean depending
% on order of AR model fit.



% Get the time series matrix from which series to be culled
[file1,path1]=uigetfile('tsm*.mat','Full-set time series matrix of chronologies:');
pf1=[path1 file1];
eval(['load ' pf1]);
if ~(exist('X')==1) | ~(exist('N')==1);
   error([pf1 'does not have X and N']);
end

% Pull year from X
yr = X(:,1);
X(:,1)=[];
[m1,n1]=size(X);

% Get longs and lats for mother set of chronologies (e.g., 669 sites)
[file5,path5]=uigetfile('crnxy?.dat','Ascii infile with full-set x,y coords');
pf5=[path5 file5];
eval(['load ' pf5]);
eval(['xycoords= ' lower(strtok(file5,'.')) ';']);
if size(xycoords,1)~=n1;
   error(['n of variables in ' pf1 ' must equal # of coords in ' pf5]);
end


% Get screen6.m output
[file2,path2]=uigetfile('*os6.mat','?p#os6.mat output file');
pf2=[path2 file2];
eval(['load ' pf2]);
if ~(exist('Lcull')==1) | ~(exist('Lneg')==1) | ~(exist('logscrn')==1);
   error([pf2 ' does not have Lcull, Lneg and logscrn']);
end



%*****************  OPTIONAL SIGN CHANGE OF SOME CHRONOLOGIES


% Note here operating on full-column matrix before culling screened series
if any(Lneg);
   ButtonName=questdlg('Sign change negative-r chronolgies?');
   switch ButtonName;
   case 'Yes';
      XN = X(:,Lneg); % submtx of chrons to be sign changed
      XN = -1.0 * XN;
      X(:,Lneg)=XN;
      ineg=1; % flag saying some chrons have been sign changed
      % Build logical pointer to culled series indicating which have been sign changed
      Ltemp=zeros(m1,1);
      Ltemp(Lcull)=(1:sum(Lcull));
      Ltemp1 = Lcull & Lneg;
      Isc = Ltemp(Ltemp1); % index pointer to sign changed series in output X
           
   case 'No';
      % No action needed
      ineg = 0;
      Isc=[];
      
   case 'Cancel';
      ineg=0;
      Isc=[];
   end
else
   ineg=0;
end
clear Ltemp Ltemp1


%**********  PROMPT FOR TIME COVERAGE OF OUTPUT MATRIX
  
% Cull rows, after prompting for time coverage
yrgo1=yr(1); % start year of full input matrix
yrsp1=yr(m1);
prompt={'Enter the start year:','Enter the end year:'};
def={yrgo1,yrsp1};
title='Desired Period for Output TSM';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
yrgo2 = str2num(answer{1});
yrsp2 = str2num(answer{2});
 


 
%***************  PULL COLUMNS FOR SELECTED SERIES

% Cull columns
X= X(:,Lcull);

% Cull long lats
xycoords=xycoords(Lcull,:);

%********************  AR MODELING AND WHITENING

ButtonName=questdlg('Do you also want a matrix of AR residuals of chronologies?');
switch ButtonName;
case 'Yes'
   arflag=1;
   %-----------  Prompt for maximum AR order to consider
   prompt={'Enter the maximum order'};
   def={2};
   title='Maximum AR order to Consider in Whitening';
   lineNo=1;
   answer=inputdlg(prompt,title,lineNo,def);
   maxord= str2num(answer{1});
   
   %------------ Fit model and whiten
   clc;
   disp('Hold on -- fitting AR models...might take a minute');
   
   [E,K1,W]=whit2(X,maxord);
   %   E (m1 x n1) the AR residuals, mean of X added back in
   %   K1 (1 x n1) order of fitted AR model
   %   W (13 x n1) information on fit:
   %	Cols 1:k1(i) - the AR coefficients
   %	Cols 7:k1(i)+6 - the corresponding standard errors
case 'No';
   arflag=0;
   E=[]; K1=[];  W=[];
case 'Cancel'
   arflag=0;
   E=[]; K1=[];  W=[];

end


%******************  CULL OUTPUT YEARS **************************

Lyr = yr>=yrgo2 & yr<=yrsp2;
yr2=(yrgo2:yrsp2)';
X=X(Lyr,:);
if arflag==1;
   E=E(Lyr,:);
end
% Note that year column still not back on X and E

nsers2 = size(X,2);


disp('Starting to build output files');

%**************** MAKE ASCII FILE OF COORDINATES

[file6,path6]=uiputfile('xy?.dat','Outfile of ascii long lats');
pf6=[path6 file6];
fid6=fopen(pf6,'w');
fmtxy='%9.3f %9.3f\n';
fprintf(fid6,fmtxy,xycoords');
fclose(fid6);


%*************** MAKE ASCII SOV AND HEADER LINE

% Original (not whitened, but possibly with a few sign-changed)
[file4,path4]=uiputfile('tvect*.dat','output SOV of original data');
pf4=[path4 file4];
fid4=fopen(pf4,'w');
string1 = sprintf('%4.0f-%4.0f, %3.0f series   ; Original (not AR Residuals)',...
   yrgo2,yrsp2,nsers2);
fprintf(fid4,'%s\n',string1);
x = X(:);
fprintf(fid4,'%6.1f\n',x);
fclose(fid4);

if arflag==1;
   % AR residuals
   [file4,path4]=uiputfile('tvect*.dat','Output SOV of AR residuals');
   pf4=[path4 file4];
   fid4=fopen(pf4,'w');
   string1 = sprintf('%4.0f-%4.0f, %3.0f series   ; AR Residuals)',...
      yrgo2,yrsp2,nsers2);
   fprintf(fid4,'%s\n',string1);
   x = E(:);
   fprintf(fid4,'%6.1f\n',x);
   fclose(fid4);
end



%*********************  MAKE .MAT OUTPUT FILE

%--------  put year col back on
X=[yr2 X];

if arflag==1;
   E=[yr2 E];
end

% Names subset char mtx
N=N(Lcull,:);

% Variables to save
saveset = ' X E N W K1 xycoords Isc';

% Save 
[file3,path3]=uiputfile('*os7.mat','?P#os7.mat output storage');
pf3=[path3 file3];
eval(['save ' pf3   saveset]); 
